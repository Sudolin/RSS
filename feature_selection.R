library(dplyr)
library(survival)
library(Boruta)
library(future.apply)
library(survivalROC)

repstress.gene<-read.csv(file = "~/Repstress/fig1/feature.selection/repstress_signature.csv")
load("~/tcga_log2TPM_dfs.Rdata")

tcgaexp_dfs <- tcgaexp_dfs %>% t() %>% 
  as.data.frame() %>%
  tibble::rownames_to_column("Sample ID") 

repstress.gene <- unique(repstress.gene$Gene)
intersect(colnames(tcgaexp_dfs),repstress.gene)

tcga_rep <- dfs %>%
  inner_join(tcgaexp_dfs[,c("Sample ID",intersect(colnames(tcgaexp_dfs),repstress.gene))],by= "Sample ID")


##  univariate cox regression--------------------------------------------------------------------
res <- data.frame()
data <- tcga_rep
genes <- colnames(data)[19:ncol(data)]

for (i in 1:length(genes)) {
  print(genes[i])
  surv =as.formula(paste('Surv(DFS, DFS_Status)~', genes[i]))
  x = coxph(surv, data = data)
  x = summary(x)
  coef=x$coef[1]
  p.value=signif(x$wald["pvalue"], digits=2)
  HR =signif(x$coef[2], digits=2);#exp(beta)
  HR.confint.lower = signif(x$conf.int[,"lower .95"], 2)
  HR.confint.upper = signif(x$conf.int[,"upper .95"],2)
  CI <- paste0("(", 
               HR.confint.lower, "-", HR.confint.upper, ")")
  res[i,1] = genes[i]
  res[i,2] = coef
  res[i,3] = HR
  res[i,4] = CI
  res[i,5] = p.value
}

names(res) <- c("gene","coef","HR","95% CI","p.value")
res.sig <- filter(res, p.value <0.01)

### Bootstrap
outTab <- NULL
surv <- tcga_rep[,c("DFS", "DFS_Status",res.sig$gene)]
rownames(surv) <- surv$`Patient ID`
for(i in 3:ncol(surv)){ # survival information (OS in this case)
  
  # customized function
  display.progress = function (index, totalN, breakN=20) {
    if ( index %% ceiling(totalN/breakN)  ==0  ) {
      cat(paste(round(index*100/totalN), "% ", sep=""))
    }
  }    
  
  display.progress(index = i, totalN = ncol(surv)) # show running progression
  gene <- colnames(surv)[i]
  Mboot <- future_replicate(1000, expr = { # bootstrap for 1,000 times
    indices <- sample(rownames(surv), size = nrow(surv) * 0.8, replace = F) # extract 80% samples at each bootsratp
    data <- surv[indices,]
    fmla1 <- as.formula(Surv(data[,"DFS"],data[,"DFS_Status"]) ~ data[,gene])
    mycox <- coxph(fmla1,data = data)
    coxResult <- summary(mycox)
    P <- coxResult$coefficients[,"Pr(>|z|)"]
  }
  )
  times <- length(Mboot[which(Mboot < 0.01)])
  outTab <- rbind(outTab,
                  cbind(gene = gene,
                        times = times))
}

outTab <- as.data.frame(outTab)

bootGene <- outTab[as.numeric(as.character(outTab$times)) > 800,] 
# Boruta ------------------------------------------------------------------
boruta =Boruta (DFS_Status ~ ., data = surv[,c("DFS_Status",share)], doTrace=2, ntree=1000, maxRuns = 1000)
table(boruta$finalDecision)
pdf("~/Repstress/fig1/feature.selection/boruta_plotImpHistory.pdf",height = 4,width = 5)
Boruta::plotImpHistory(boruta)
dev.off()
plot(boruta)