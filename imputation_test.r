#function to test significance across columns/different variables for IMPUTATION PANELS
imputation_wilcox_test <- function(x,class_factor,score, alternative = c("two.sided", "less", "greater")) {
    #panels is an array of selected panels to test
    #score is what we want to test: INFO, R2, CONCORDANCE etc..
    #class_factor is an array of column names, used to subset X
    alternative <- match.arg(alternative)
    f1 <- class_factor[1]
    f2 <- class_factor[2]
    f1_val <- sort(unique(x[,f1]))
    f1_ans_all <- NULL
    for (i in f1_val){
        x_f1 <- x[x[,f1] %in% i,]
        ans <- pairwise.wilcox.test(x_f1[,score],x_f1[,f2],alternative=alternative)
        f1_ans <- as.data.frame(ans$p.value)
        f1_ans[,f1] <- i
        f1_ans$alternative <- alternative
        f1_ans_all <- rbind(f1_ans_all,f1_ans)
    }
    f1_ans_all
}