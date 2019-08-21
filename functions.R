
prejoin <- function(dt, ext=250, width=NULL){
  if(is.null(width)){ # extend CI ends by fixed amount
    dt[,CI_start_ext := CI_start -ext]
    dt[,CI_stop_ext := CI_stop +ext]

  }else{ # set total peak width as set value around center
    dt[,CI_start_ext := center_start - width/2]
    dt[,CI_stop_ext := center_stop + width/2]
  }
  setkey(dt, chr, CI_start_ext, CI_stop_ext)
  return(dt)
}

trim_peaks <- function(dt, max_len=1000){
  dt[,width := CI_stop-CI_start]
  dt[,CI_start_ext := as.integer(CI_start_ext)]
  dt[,CI_stop_ext := as.integer(CI_stop_ext)]
  dt[,CI_start_ext := as.integer(CI_start + max(0,(width/2)-(max_len/2))), by=c("chr", "CI_start", "CI_stop")]
  dt[,CI_stop_ext := as.integer(CI_stop - max(0,(width/2)-(max_len/2))), by=c("chr", "CI_start", "CI_stop")]
  setkey(dt, chr, CI_start_ext, CI_stop_ext)
  return(dt)
}


fraction_overlapN_old <- function(dt,x,y,n=100){
  setorder(dt, enrichment)
  binning_vector <- cut(1:nrow(dt), n)
  plot(1:n/n, tapply(!is.na(dt[,get(y)]), binning_vector, mean),
       xlab=paste("Fraction of",x,"Peaks  (ordered low to high enrichment)"),
       ylab=paste("Fraction of",x,"peaks that overlap",y,"peaks"))
}



fraction_overlapN <- function(dt,x,y, n=25, filters=c("All"), returndf=F){
  setorder(dt, enrichment)
  if(is.integer(dt$enrichment)){
    dt$enrichment_bin <- dt$enrichment
    warning("Using enrichment as bin (integer type detected)")
  }else{
    dt$enrichment_bin <- cut(1:nrow(dt), n)
  }


  All=TRUE

  tmp <- dt[eval(parse(text=filters[1])),.(mean_enrichment=mean(enrichment),
                                           num=.N,
                                           Overlap=mean(!is.na(get(y))),
                                           Subset=filters[1],
                                           #sem=sd(!is.na(get(y)))/sqrt(length(get(y))),
                                           sep=sqrt((mean(!is.na(get(y))) * (1- mean(!is.na(get(y)))))/length(get(y)))),by=enrichment_bin]

  #dt[,mean(!is.na(get(y)))]

  if(returndf){
    return(tmp)
  }else{
    return(
      ggplot(tmp, aes(mean_enrichment, Overlap, colour=Subset, group=Subset)) +
        #geom_point(size=2) +
        geom_hline(yintercept = dt[,mean(!is.na(get(y)))], colour="blue", linetype="dashed") +
        geom_pointrange(aes(ymax = Overlap + 2*sep, ymin = Overlap - 2*sep)) +
        theme_minimal() +
        scale_color_brewer(palette = 'Set1') +
        xlab(paste("Mean",x,"Enrichment")) +
        scale_x_log10() +
        ylab(paste("Fraction of ",x,"peaks\n that overlap",y,"peaks"))
    )}
}


randomise_positions <- function(base, shift){
  base_chr_lengths <- base[,.("chr_len"=max(CI_stop_ext)),by=chr]
  setkey(base_chr_lengths, chr)

  base_shifted <- base[base_chr_lengths, on="chr"]

  base_shifted[,max_right_shift := min(chr_len,CI_stop_ext+shift) - CI_stop_ext, by=c("chr", "CI_start_ext", "CI_stop_ext")]
  base_shifted[,max_left_shift := CI_start_ext - max(0,CI_start_ext-shift), by=c("chr", "CI_start_ext", "CI_stop_ext")]

  set.seed(42)
  #shifts <- sample(-shift:shift, nrow(base_shifted), replace = T)
  base_shifted[,shift := sample(-max_left_shift:max_right_shift,1), by=c("chr", "CI_start_ext", "CI_stop_ext")]

  base_shifted[,CI_start_ext := CI_start_ext + shift]
  base_shifted[,CI_stop_ext := CI_stop_ext + shift]
  setkey(base_shifted, chr, CI_start_ext, CI_stop_ext)

  return(base_shifted)
}

fraction_overlapN2 <- function(x,y, x1, y1, n=25, shift=1e5, filters=c("All"), returndf=F, corrected=F){

  base <- x[,.(chr, CI_start_ext, CI_stop_ext, enrichment)]
  test <- y[,.(chr, CI_start_ext, CI_stop_ext)] # enrichment

  setkey(base, chr, CI_start_ext, CI_stop_ext)
  setkey(test, chr, CI_start_ext, CI_stop_ext)

  base_shifted <- randomise_positions(base, shift)

  base$test_enrichment <- foverlaps(base, test, mult="first")$CI_start_ext # enrichment
  base_shifted$test_enrichment <- foverlaps(base_shifted, test, mult="first")$CI_start_ext # enrichment

  a <- fraction_overlapN(base, "base", "test_enrichment", n=n, returndf = T)
  b <- fraction_overlapN(base_shifted, "base_shifted", "test_enrichment",n=n, returndf = T)

  a$Type <- "Observed"
  b$Type <- "Randomised"

  tmp <- rbind(a,b)

  if(corrected){
    tmp <- dcast(tmp, enrichment_bin + mean_enrichment + num + Subset ~ Type, value.var = "Overlap")
    tmp[,Overlap := (1 - ( (1-Normal) / (1-Randomised) ))]
    tmp[,Type := "Corrected"]
    tmp[,sep := sqrt((Overlap * (1- Overlap))/num), by=enrichment_bin]
  }


  if(returndf){
    return(tmp)
  }else{

    p <- ggplot(tmp, aes(mean_enrichment, Overlap, colour=Type, group=Type))

    if(!corrected){
      p <- p + geom_hline(yintercept = base[,mean(!is.na(test_enrichment))], colour="#E41A1C", linetype="dashed") +
        geom_hline(yintercept = base_shifted[,mean(!is.na(test_enrichment))], colour="#377EB8", linetype="dashed")
    }

    return(p +
             geom_pointrange(aes(ymax = Overlap + 2*sep, ymin = Overlap - 2*sep)) +
             theme_minimal() +
             scale_color_brewer(palette = 'Set1') +
             xlab(paste("Mean",x1,"Enrichment")) +
             scale_x_log10() +
             ylim(-0.01,1.01) +
             ylab(paste("Fraction of ",x1,"peaks\n that overlap",y1,"peaks"))
    )}
}



fraction_overlapP <- function(dt,x,y,n=100){
  setorder(dt, -pvalue)
  binning_vector <- cut(1:nrow(dt), n)
  plot(1:n/n, tapply(!is.na(dt[,get(y)]), binning_vector, mean),
       xlab=paste("Fraction of",x,"Peaks  (ordered large to small pvalue)"),
       ylab=paste("Fraction of",x,"peaks that overlap",y,"peaks"))
}


loadBed <- function(file, ext=250){
  Trim28 <- fread(file)
  names(Trim28)[1:3] <- c("chr","CI_start","CI_stop")
  Trim28 <- prejoin(Trim28, ext)
  return(Trim28)
}

overlapGam <- function(x, y, shift=NULL, range=c(0.01,50), returngam=F){

  base <- x[,.(chr, CI_start_ext, CI_stop_ext, enrichment)]
  test <- y[,.(chr, CI_start_ext, CI_stop_ext)] # enrichment

  setkey(base, chr, CI_start_ext, CI_stop_ext)
  setkey(test, chr, CI_start_ext, CI_stop_ext)


  if(!is.null(shift)){
    base <- randomise_positions(base, shift)
  }

  base$test_enrichment <- foverlaps(base, test, mult="first")$CI_start_ext # enrichment

  base$overlap <- !is.na(base$test_enrichment)

  gam_fit <- gam(data=base, overlap ~ s(enrichment), family="quasibinomial")

  if(returngam){
    return(gam_fit)
  }else{
    n=500
    x_new <- seq(range[1], range[2], length.out = n)
    gam_pred <- data.table(enrichment = x_new,
                           data.frame(predict(gam_fit, data.frame(enrichment = x_new), se.fit=TRUE, type="response")))
    gam_pred <- cbind(gam_pred, Type="GAM")
    setnames(gam_pred, c("mean_enrichment","Overlap","se","Type"))
    return(gam_pred)
  }

}
