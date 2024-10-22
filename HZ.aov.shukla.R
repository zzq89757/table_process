###updated by gqyao 2017.9.27
###Modified for treat missing data. replace missing data with y= u + v + rep for each site.
HZ.aov.shukla <-function(filename){
    sub_fun_path <- "/mnt/ntc_data/wayne/Repositories/table_process/"
    sub_fun_path <- paste(sub_fun_path, "shareFun4.1.R", sep="/")
    source(sub_fun_path)
    library(reshape)
    avo.file.path<-paste("/mnt/ntc_data/wayne/Repositories/table_process/",filename,".txt",sep="")
    yield <- read.delim(avo.file.path, header=T, sep="\t",fileEncoding="UTF-8-BOM", stringsAsFactors=F )
    yield <- ZeroToNA(yield)
    ## set the Mean Square Error, and df from the Anov results.
    aov.result <-hz.aov_tot(yield)[[1]]
    n.mse <- aov.result[[1]][5,3]
    n.mse_df <- aov.result[[1]][5,1]
    ## covert number of locations and variety to environments and genotypes respectively.
    num.env <- length(unique(yield$location))
    num.gen <- length(unique(yield$variety))
    ## number of replicates for each variety in a given environment
    num.rep <- ncol(yield)-2
    ## estimate missingdata. Coded by Guoqi 
    if(nMissdata(yield) > 0) yield <- lstSqEstMissdata(yield)
    ##calculate Shukla's stability variance when missing data estimated 
    yield.mat <- yield[,-c(1:2)]
    row.means <- rowMeans(yield.mat)
    df.avg.yield <- data.frame(location =yield[,1], variety=yield[,2], avg.yield=row.means)
    avg.yield <- cast(df.avg.yield, variety~location, value.var="avg.yield")
    print(avg.yield)
    rownames(avg.yield ) <- avg.yield [,1]
    avg.yield <- avg.yield[,-1]
    #calculate Shukla's stability variance, recoded by gqyao
    v.loc_mean <- colMeans(avg.yield)
    v.var_mean <- rowMeans(avg.yield)
    y_ave <- mean(v.var_mean)
    # caculate shukla parameter
    mt.loc_ave <- v.loc_mean
    for (i in c(1:(num.gen-1))) {mt.loc_ave = rbind(mt.loc_ave, v.loc_mean)}
    mt.uij <- avg.yield - mt.loc_ave
    v.ui_ave <- rowMeans(mt.uij)
    mt.ui_ave <- v.ui_ave
    print(mt.uij)
    for(i in c(1:(num.env-1))){mt.ui_ave <- cbind(mt.ui_ave, v.ui_ave)}
    mt.uij.Minusui_ave.sq <- (mt.uij - mt.ui_ave)^2
    v.uij.Minusui_ave.sq.rsum <- rowSums(mt.uij.Minusui_ave.sq)
    v.uij.Minusui_ave.sq.sum <- sum(v.uij.Minusui_ave.sq.rsum)
    v.var_name <- rownames(avg.yield)
    v.shukla_var <- (num.gen*(num.gen-1)*v.uij.Minusui_ave.sq.rsum - v.uij.Minusui_ave.sq.sum)/(num.env-1)/(num.gen-1)/(num.gen -2)
    v.shukla_var[v.shukla_var < 0] <- 0 #check v.shukla_var < 0
    v.F <- v.shukla_var/(n.mse/num.rep)
    v.P = pf(v.F, num.env-1, n.mse_df, lower.tail = F)
    v.interact_var <-  v.shukla_var - n.mse/num.rep
    v.interact_var[v.interact_var < 0] <- 0 #check  v.interact_var < 0
    if (sum((v.var_mean <= 0) * 1) > 0 ) stop("a variety mean yield <= 0")
    v.CV <- sqrt(v.shukla_var) *100/v.var_mean
    output <- data.frame(variety=v.var_name, df=num.env-1, Shukla.Var=v.shukla_var, F=v.F, 
			  P=v.P, interaction.var= v.interact_var, variety.mean=v.var_mean,
			 Shukla.CV=v.CV, stability.rank=order(v.CV))
    ##Bartlett test
    n.Sp_sq <- sum((num.env-1)*v.shukla_var)/(num.gen*num.env-num.gen)
    n.chi_sq <- ((num.gen*num.env-num.gen)*log(n.Sp_sq) - sum((num.env-1)*log(v.shukla_var)))/(1 + (1/(num.env-1) * num.gen - 1/(num.env*num.gen - num.gen))/3/(num.gen-1))
    n.p.Bartlett <- pchisq(n.chi_sq, num.gen-1,lower.tail = F)
    ##Shukla variance F test
    v.shukla_var_ordered <- v.shukla_var[order(v.shukla_var,decreasing = T)]
    mt.shukla_var_F <- v.shukla_var_ordered %*% t(1/v.shukla_var_ordered)
    rownames(mt.shukla_var_F) <- colnames(mt.shukla_var_F)
    mt.shukla_var_F[lower.tri(mt.shukla_var_F, diag = T)] <- NA
    mt.shukla_var_F_P <- pf(mt.shukla_var_F,num.env-1,num.env-1,lower.tail=F)
    print(mt.ui_ave)
    v.05_comp <- create_FComp_vec(0.05, F, mt.shukla_var_F_P)
    v.01_comp <- create_FComp_vec(0.01, T, mt.shukla_var_F_P)
    output_F_comp <- data.frame(variety =  rownames(mt.shukla_var_F_P), F05_level = v.05_comp, F01_level = v.01_comp)
    output_merge <- merge(output,output_F_comp,by = c("variety"),all = T)
    write.table(output_merge, "Analysis of stability by Shukl method.txt", sep ="\t", quote = FALSE, row.names = FALSE)
    return(list(output_merge, n.p.Bartlett))
}

###########################################################################################################
### A function to creat and return comparison vector 
create_FComp_vec <- function(n.threshold, b.cap = F, mt.comp){
#n.threshold: the standard value (including)
#b.cap: T to use capital letters, F to use lowercase letters
#mt.comp: up tridialog intact comparison information with row to be compared with column
    num.gen <- ncol(mt.comp)
    v.comp = c()
    n.grp = 1
    i = 1
    while (i < num.gen){
	if (is.null(v.comp[i])){
	    if (n.grp > 26){# when group no is larger than 26 ,the max number of letter table.
		if (b.cap){
		    v.comp[i]= paste(LETTERS[n.grp%%26],as.character(n.grp%/%26), sep = "")
		}else{
		    v.comp[i]= paste(letters[n.grp%%26],as.character(n.grp%/%26), sep = "")
		}
	    }else{
		if (b.cap){
		    v.comp[i] = LETTERS[n.grp]
		}else{
		    v.comp[i] = letters[n.grp]
		}
	    }
	}
	for(j in ((i+1):num.gen)){
        # print(i)
	    if (is.nan(mt.comp[i,j])) stop("Can not compare to 0/0!")
	    if (mt.comp[i,j] > n.threshold){
		v.comp[j] = v.comp[i]
	    }
	    else{
		n.grp = n.grp + 1
		if (n.grp > 26){
		    if (b.cap){
			v.comp[j] = paste(LETTERS[n.grp%%26],as.character(n.grp%/%26), sep = "")
		    }else{
			v.comp[j] = paste(letters[n.grp%%26],as.character(n.grp%/%26), sep = "")
		    }
		}else{
		    if (b.cap){
			v.comp[j] = LETTERS[n.grp]
		    }else{
			v.comp[j] = letters[n.grp]
		    }
		}
		#lookup backup
		for(k in ((j -1):1)){
		    if (is.nan(mt.comp[k,j])) stop("Can not compare to 0/0!")
		    if (mt.comp[k,j] <= n.threshold) break
		    v.comp[k] = paste(v.comp[k], v.comp[j], sep= " ")
		}
		#lookupforword
		i = j
		break
	    }
	}
	if (j == num.gen) break
    }
    return (v.comp)
}

HZ.aov.shukla("AOV_breed_120241021145105772")

####END######
