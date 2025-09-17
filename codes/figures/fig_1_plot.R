cd_path <- "D:/ITB_S-1/Yr_4/KP/connection_analysis"
setwd(cd_path)

df_exp <- as.data.frame(
  read.csv(
    "expressionAnalysis_merge.csv",sep=',',header=TRUE, 
    fill=TRUE,quote=""
  )
)
clinical_exp <- round(df_exp[2:40]*1000)
laboratory_exp <- round(df_exp[41:57]*1000)

gene_ids <- as.vector(unlist(df_exp[1]))
len_df <- length(gene_ids)

clinical_readsum <- as.vector(colSums(clinical_exp))
laboratory_readsum <- as.vector(colSums(laboratory_exp))

min_samp1 <- round(median(clinical_readsum*0.8))
min_samp2 <- round(median(laboratory_readsum*0.8))
min_samp <- min(c(min_samp1, min_samp2))
clinical_exp <- clinical_exp[,which(clinical_readsum>=min_samp)]
laboratory_exp <- laboratory_exp[,which(laboratory_readsum>=min_samp)]

#min_samp1 <- quantile(clinical_filter,probs=0.25) - 1.5*(quantile(clinical_filter,probs=0.75)-quantile(clinical_filter,probs=0.25))
#min_samp2 <- quantile(laboratory_filter,probs=0.25) - 1.5*(quantile(laboratory_filter,probs=0.75)-quantile(laboratory_filter,probs=0.25))
#min_samp <- min(c(min_samp1, min_samp2))
#clinical_expression <- clinical_expression[,which(clinical_filter>=min_samp)]
#laboratory_expression <- laboratory_expression[,which(laboratory_filter>=min_samp)]

sub_samp1 <- round(median(clinical_readsum*0.8))
sub_samp2 <- round(median(laboratory_readsum*0.8))
sub_samp <- min(c(sub_samp1,sub_samp2))
for (x in 1:ncol(clinical_exp)) {
  clinical_exp[,x] <- as.vector(rrarefy(clinical_exp[,x], sample=50000))
}
for (x in 1:ncol(laboratory_exp)) {
  laboratory_exp[,x] <- as.vector(rrarefy(laboratory_exp[,x], sample=50000))
}

p_values <- vector(mode="numeric",length=len_df)
for (i in 1:len_df) {
  tryCatch(
    expr = {
      clinical_group <- as.numeric(clinical_exp[i,])
      laboratory_group <- as.numeric(laboratory_exp[i,])
      MW_result <- wilcox.test(
        x=clinical_group,y=laboratory_group,exact=FALSE,
        paired=FALSE,alternative="two.sided",conf.level=0.95
      )
      p_values[i] <-  MW_result$p.value
    },
    error = function(e) {
      p_values[i] <-  NA
    }
  )
}
p_values <- p.adjust(p_values,method='BH')
log_p <-  -log10(p_values)


func_mean <- function(data,indices) {
  data_samp <-data[indices]
  return(mean(data_samp))
}

clinical_means <- vector(mode='numeric',length=len_df)
for (i in 1:len_df) {
  bootstrap_result <- boot(as.vector(unlist(clinical_exp[i,])),statistic=func_mean,R=100)
  clinical_means[i] <- bootstrap_result$t0
}

laboratory_means <- vector(mode='numeric',length=len_df)
for (i in 1:len_df) {
  bootstrap_result <- boot(as.vector(unlist(laboratory_exp[i,])),statistic=func_mean,R=100)
  laboratory_means[i] <- bootstrap_result$t0
}

log2_fold <- log2(clinical_means/laboratory_means)

#remove inf -inf nan
keep <- is.finite(log_p) & is.finite(log2_fold)
log_p <- log_p[keep]
log2_fold <- log2_fold[keep]
gene_ids <- gene_ids[keep]

len_df <- length(gene_ids)

p_threshold =-log10(0.05)
fold_threshold = 2.5

pf <- abs(log_p)>p_threshold & abs(log2_fold)>fold_threshold
p <- abs(log_p)>p_threshold & abs(log2_fold)<=fold_threshold
f <- p_threshold>=abs(log_p) & abs(log2_fold)>fold_threshold
ns <- p_threshold>=abs(log_p) & fold_threshold>=abs(log2_fold)
pf_high <- log_p>p_threshold & log2_fold>fold_threshold


is_high <- vector(mode='logical',length=len_df)
is_high[pf_high] <- TRUE

category <- vector(mode='character',length=length(gene_ids))
category[pf] <- "pf"
category[f] <- "f"
category[p] <- "p"
category[ns] <- "ns"

df <- data.frame(gene_id=gene_ids,log_p=log_p,log2_fold=log2_fold,category=category,is_high=is_high)

gene_ids_high <- vector(mode='character',length=len_df)
gene_ids_high[is_high] <- gene_ids[is_high]

ggplot(df,aes(x=log_p,y=log2_fold,alpha=1,color=category)) + 
  geom_point(show.legend=TRUE) +
  scale_color_manual(
    values = c(
      "pf"= "red",
      "f"  = "green",
      "p" = "blue",
      "ns" = "black"
    )
  ) +
  geom_text(size=2,aes(label=gene_ids_high)) + 
  geom_hline(yintercept=c(-fold_threshold,fold_threshold)) +
  geom_vline(xintercept=-log10(0.05))

high_only <- as.data.frame(gene_ids[which(is_high)])
write.csv(high_only,"high_genes.csv")

write.csv(df,"differentialAnalysis.csv")
