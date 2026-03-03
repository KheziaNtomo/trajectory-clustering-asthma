# ==============================================================================
# functions.R — Utility Functions
# ==============================================================================
# Description: Shared helper functions for plotting (PCA, biplots, ROC curves),
#              enrichment analysis, and data formatting.
#
# Functions:
#   - CreateScorePlot()      PCA score plot
#   - CreateBiPlot()         PCA biplot (base R)
#   - CreateBiPlot2()        PCA biplot (ggplot2)
#   - SaveExcelWithSuperscripts()  Excel export with formatted superscripts
#   - rocfoldslogistic()     Cross-validated ROC curves for logistic regression
#   - PlotROCCurves()        Plot ROC curves with confidence bands
#   - perform_enrichment()   Gene set enrichment (hypergeometric/Fisher)
#   - clean_enrichment()     Post-process enrichment results
#   - geneid_to_protein()    Gene ID to protein name mapping
# ==============================================================================

var_orders = list(c("Nasal polyps (Y/N)","HDM-sensitive (Y/N)","FVC","FEV1","Neutrophil", "Eosinophil","Case (Y/N)"),
                  c("Severe exacerbation (Y/N)", "Exacerbation (Y/N)", "ACQ5",  "Age of onset","Current OCS use (Y/N)","Severe (Y/N)"))

protein_levels = c("Periostin", "CCL-18", "IL-18","CCL-5","ITAC","TARC","IL-10",
                   "IL-16","MIP-1a","PAPP-A","IL-17","TSLP","SP-A","CTACK","KL-6","MPO", 
                   "EDN","IL-6", "TNFa","MIG","CCL-20","IP-10","IL-4","IL-5")

names_dict<-list(Nasalpolyps = "Nasal polyps (Y/N)",  
                 HDMsensitive = "HDM-sensitive (Y/N)", 
                 FVC = "FVC",          
                 FEV1 = "FEV1",       
                 Neutrophil = "Neutrophil",   
                 Eosinophil = "Eosinophil",   
                 case_control = "Case (Y/N)",
                 Severeexacerbation = "Severe exacerbation (Y/N)", 
                 Exacerbation = "Exacerbation (Y/N)",       
                 ACQ5 = "ACQ5",               
                 Ageofonset = "Age of onset",         
                 OCS = "Current OCS use (Y/N)",                
                 Severity =  "Severe (Y/N)")


# Randomly sample from truncated Gaussian distribution
rtruncnorm = function(n, mean = 0, sd = 1, min = -Inf, max = Inf){
  set.seed(202204)
  if (min > max) stop('Error: Truncation range is empty')
  x = runif(n, pnorm(min, mean, sd), pnorm(max, mean, sd))
  qnorm(x, mean, sd)
}

# Create PCA score plot
CreateScorePlot=function(mypca, filename=NULL, type, type.colours, comp=NULL, pch=19,
                         ellipse = FALSE, legend = TRUE, wide = FALSE){
  if (is.null(comp)){
    comp=matrix(c(1,2,1,3,2,3), byrow=TRUE, ncol=2)
  }
  if (legend){
    extra = ceiling(length(type.colours)/25)
  } else {
    extra = 0
  }
  if (!is.null(filename)){
    pdf(paste0(filename), width=15+1.7*extra, height=5) 
  }
  if (legend){
    par(mar = c(5,5,1,1))
    layout(matrix(c(1,2,3,4), 1, 4, byrow = TRUE), widths=c(3,3,3,extra))
  } else {
    par(mar = c(5,5,1,1), mfrow = c(1,3))
  }
  ev=mypca$eig[,2]
  for (k in 1:nrow(comp)){
    xcomp=comp[k,1]
    ycomp=comp[k,2]
    S = mypca$ind$coord[,c(xcomp,ycomp)]
    plot(S, pch=pch, cex=2, las=1, type = "n",
         col=type.colours[type],
         xlab=paste0("Comp ",xcomp," (", round(ev[xcomp], digits=2), "% e.v.)"),
         ylab=paste0("Comp ",ycomp," (", round(ev[ycomp], digits=2), "% e.v.)"),
         cex.lab = 1.5)
    points(S, pch=pch, cex=1.2, col=type.colours[type])
    families = unique(type)
    if (ellipse){
      for (f in 1:length(families)){
        tmpmat=S[type==families[f],]
        if (!is.null(nrow(tmpmat))){
          if (nrow(tmpmat)>2){
            lines(ellipse::ellipse(cov(tmpmat), centre = c(mean(tmpmat[,1]),mean(tmpmat[,2])), level = 0.95),
                  lty = 3, col = type.colours[unique(type)[f]])
          }
        }
        # else {
        #   plotrix::draw.ellipse(tmpmat[1], tmpmat[2], 0.5, 0.5,lty = 3, border = type.colours[f])
        # }
      } 
    }
    abline(v=axTicks(1), lty=3, col="grey")
    abline(h=axTicks(2), lty=3, col="grey")
  }
  if (legend){
    par(mar = c(0,1,1,0))
    plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n', xlab = "", ylab = "")
    legend("left", col=type.colours, ncol = ceiling(length(type.colours)/25),
           pch=pch, pt.cex=2, legend=names(type.colours), bty = "n", cex = 1.5)
  }
  if (!is.null(filename)){
    print("Saved to filename")
    dev.off()
  }
}


CreateBiPlot2 <- function(mypca, filename = NULL, type, type.colours, var.colours, 
                         comp = NULL, pch = 19, legend = TRUE, magnification = 1) {
  if (is.null(comp)) {
    comp <- matrix(c(1, 2, 1, 3, 2, 3), byrow = TRUE, ncol = 2)
  }
  ev <- mypca$eig[, 2]
  plots <- list()
  
  for (k in 1:nrow(comp)) {
    xcomp <- comp[k, 1]
    ycomp <- comp[k, 2]
    S <- mypca$ind$coord[, c(xcomp, ycomp)] / magnification
    loadings <- sweep(mypca$var$coord, 2, sqrt(mypca$eig[1:ncol(mypca$var$coord), 1]), FUN = "/")
    
    plot_data <- data.frame(S, type = type)
    colnames(plot_data) <- c("x", "y", "type")
    loadings_data <- data.frame(loadings[, c(xcomp, ycomp)])
    colnames(loadings_data) <- c("x", "y")
    loadings_data$labels <- rownames(loadings)
    
    p <- ggplot() +
      geom_point(data = plot_data, aes(x = x, y = y, color = type), shape = pch, size = 2, alpha = 0.25) +
      geom_segment(data = loadings_data, aes(x = 0, y = 0, xend = x, yend = y), 
                   arrow = arrow(length = unit(0.1, "inches")), color = var.colours) +
      geom_text(data = loadings_data, aes(x = x + sign(x) * 0.01, y = y + sign(y) * 0.01, label = labels),
                color = darken(var.colours, 0.25), hjust = 1, vjust = 1) +
      scale_color_manual(values = type.colours) +
      labs(x = paste0("Comp ", xcomp, " (", round(ev[xcomp], digits = 2), "% e.v.)"),
           y = paste0("Comp ", ycomp, " (", round(ev[ycomp], digits = 2), "% e.v.)")) +
      theme_minimal(base_size = 15) +
      theme(legend.position = "none")
    
    plots[[k]] <- p
  }
  
  if (legend) {
    legend_plot <- ggplot() +
      geom_point(aes(x = 1, y = 1, color = factor(names(type.colours))), shape = pch, size = 5) +
      scale_color_manual(values = type.colours) +
      theme_void() +
      theme(legend.position = "left", legend.text = element_text(size = 15)) +
      guides(color = guide_legend(override.aes = list(size = 5), ncol = ceiling(length(type.colours) / 25)))
    
    plots <- c(plots, list(legend_plot))
  }
  
  if (!is.null(filename)) {
    pdf(paste0(filename), width = 15 + 1.7 * ceiling(length(type.colours) / 25), height = 5)
    grid.arrange(grobs = plots, ncol = 4)
    dev.off()
  }
  
  return(plots)
}


CreateBiPlot=function(mypca, filename=NULL, type, type.colours, var.colours,
                      comp=NULL, pch=19, legend = TRUE, magnification = 1){
  if (is.null(comp)){
    comp=matrix(c(1,2,1,3,2,3), byrow=TRUE, ncol=2)
  }
  if (legend){
    extra = ceiling(length(type.colours)/25)
  } else {
    extra = 0
  }
  if (!is.null(filename)){
    pdf(paste0(filename), width=15+1.7*extra, height=5) 
  }
  if (legend){
    par(mar = c(5,5,1,1))
    layout(matrix(c(1,2,3,4), 1, 4, byrow = TRUE), widths=c(3,3,3,extra))
  } else {
    par(mar = c(5,5,1,1), mfrow = c(1,3))
  }
  ev=mypca$eig[,2]
  for (k in 1:nrow(comp)){
    xcomp=comp[k,1]
    ycomp=comp[k,2]
    S = mypca$ind$coord[,c(xcomp,ycomp)]/magnification
    loadings=sweep(mypca$var$coord,2,sqrt(mypca$eig[1:ncol(mypca$var$coord),1]),FUN="/")
    plot(S, pch=pch, cex=2, las=1, type = "n",
         col=type.colours[type],
         ylim = range(c(S[,2],loadings[,ycomp]+sign(loadings[,ycomp])*0.1)),
         xlim = range(c(S[,1],loadings[,xcomp]+sign(loadings[,xcomp])*0.1)),
         xlab=paste0("Comp ",xcomp," (", round(ev[xcomp], digits=2), "% e.v.)"),
         ylab=paste0("Comp ",ycomp," (", round(ev[ycomp], digits=2), "% e.v.)"),
         cex.lab = 1.5)
    points(S, pch=pch, cex=1.2, col=alpha(type.colours[type],0.25))
    points(loadings[,c(xcomp,ycomp)],cex=0.1, pch=19,col=var.colours)
    arrows(x0=rep(0, nrow(loadings)), y0=rep(0, nrow(loadings)),
           x1=loadings[,xcomp], y1=loadings[,ycomp], length=0.1, col=var.colours, lwd=1.5)
    abline(v=axTicks(1), lty=3, col="grey")
    abline(h=axTicks(2), lty=3, col="grey")
    text(loadings[,xcomp]+sign(loadings[,xcomp])*0.01, loadings[,ycomp]+sign(loadings[,ycomp])*0.01,
         labels=rownames(loadings), cex=1,
         col=darken(var.colours,0.25))
  }
  if (legend){
    par(mar = c(0,1,1,0))
    plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n', xlab = "", ylab = "")
    legend("left", col=type.colours, ncol = ceiling(length(type.colours)/25),
           pch=pch, pt.cex=2, legend=names(type.colours), bty = "n", cex = 1.5)
  }
  if (!is.null(filename)){
    print("Saved to filename")
    dev.off()
  }
}

SaveExcelWithSuperscripts=function(dt, filename){
  wb <- openxlsx::createWorkbook() # create workbook
  openxlsx::addWorksheet(wb, sheetName = "data") # add sheet
  openxlsx::writeData(wb, sheet=1, x=dt, xy=c(1, 1)) # write data on workbook
  
  for(i in grep("\\_\\[([A-z0-9\\s\\-]*)\\]", wb$sharedStrings)){
    # if empty string in superscript notation, then just remove the superscript notation
    if(grepl("\\_\\[\\]", wb$sharedStrings[[i]])){
      wb$sharedStrings[[i]] <- gsub("\\_\\[\\]", "", wb$sharedStrings[[i]])
      next # skip to next iteration
    }
    
    # insert additioanl formating in shared string
    wb$sharedStrings[[i]] <- gsub("<si>", "<si><r>", gsub("</si>", "</r></si>", wb$sharedStrings[[i]]))
    
    # find the "_[...]" pattern, remove brackets and udnerline and enclose the text with superscript format
    wb$sharedStrings[[i]] <- gsub("\\_\\[([A-z0-9\\s\\-]*)\\]", "</t></r><r><rPr><vertAlign val=\"superscript\"/></rPr><t xml:space=\"preserve\">\\1</t></r><r><t xml:space=\"preserve\">", wb$sharedStrings[[i]])
  }
  openxlsx::saveWorkbook(wb, file=filename, overwrite = TRUE)
}

rocfoldslogistic=function(xdata, ydata, M=5, niter=1){
  t0=Sys.time()
  sensitivity=specificity=AUC=NULL
  for (i in 1:niter){
    # scv=sample(rep(1:M,length.out=nrow(xdata)), size=nrow(xdata))
    tmp=split(1:length(ydata), f=as.factor(ydata))
    cacofolds=lapply(tmp, FUN=function(x){sample(rep(1:M,length.out=length(x)), size=length(x))})
    scv=rep(NA,length(ydata))
    for (k in 1:2){
      scv[tmp[[k]]]=cacofolds[[k]]
    } # balanced folds
    beta=NULL
    fitted=NULL
    for (k in 1:M){ # CV
      xdata_train=xdata[scv!=k,,drop=FALSE]
      y_train=ydata[scv!=k]
      xdata_test=xdata[scv==k,,drop=FALSE]
      colnames(xdata_train)=colnames(xdata_test)=colnames(xdata)
      y_test=ydata[scv==k]
      f=paste('y_train ~', paste(colnames(xdata), collapse = ' + '))
      # print(f)
      model=glm(as.formula(f), data = as.data.frame(xdata_train), family = binomial)
      beta=c(beta, coefficients(model)[-1])
      fitted=c(fitted, predict(model, as.data.frame(xdata_test)))
    }
    roc_iter=roc(response=ydata, predictor=fitted[names(ydata)])
    sensitivity=cbind(sensitivity, roc_iter$sensitivities)
    specificity=cbind(specificity, roc_iter$specificities)
    AUC=c(AUC, roc_iter$auc[1])
  }
  t1=Sys.time()
  print(t1-t0)
  return(list(sensitivity=sensitivity, specificity=specificity, AUC=AUC, 
              nobs=nrow(xdata), ncases=sum(as.character(ydata)=="1")))
}


PlotROCCurves=function(mymodels, niter=100){
  model=mymodels[[1]]
  #model=eval(parse(text=mod))
  plot(1-apply(model$specificity, 1, mean), apply(model$sensitivity, 1, mean), type='l', lwd=2, col=colors[1],
       xlab='False Positive Rate', ylab='True Positive Rate', las=1, cex.lab=1.3,
       panel.first=abline(0,1,lty=3))
  
  for (i in 1:length(mymodels)){ # area between most extreme points (quantiles)
    model=mymodels[[i]]
    #model=eval(parse(text=mod))
    polygon(c(1-apply(model$specificity, 1, FUN=function(x){sort(x)[quantiles[1]*niter]}), 
              rev(1-apply(model$specificity, 1, FUN=function(x){sort(x)[quantiles[2]*niter]}))),
            c(apply(model$sensitivity, 1, FUN=function(x){sort(x)[quantiles[1]*niter]}), 
              rev(apply(model$sensitivity, 1, FUN=function(x){sort(x)[quantiles[2]*niter]}))),
            col=tcolors[i], border=NA)
  }
  
  for (i in 1:length(mymodels)){ # curve with most extreme points
    model=mymodels[[i]]
    #model=eval(parse(text=mod))
    qs=sort.list(model$AUC)[c(quantiles[1]*niter, quantiles[2]*niter)]
    lines(1-apply(model$specificity, 1, FUN=function(x){sort(x)[quantiles[1]*niter]}),
          apply(model$sensitivity, 1, FUN=function(x){sort(x)[quantiles[1]*niter]}), type='l', lwd=0.3, col=colors[i])
    lines(1-apply(model$specificity, 1, FUN=function(x){sort(x)[quantiles[2]*niter]}),
          apply(model$sensitivity, 1, FUN=function(x){sort(x)[quantiles[2]*niter]}), type='l', lwd=0.3, col=colors[i])
  }
  
  for (i in 1:length(mymodels)){
    model=mymodels[[i]]
    #model=eval(parse(text=mod))
    lines(1-apply(model$specificity, 1, mean), apply(model$sensitivity, 1, mean), type='l', lwd=2, col=colors[i])
  }
}


perform_enrichment <- function(genes_of_interest, pathways, test_type = "hypergeomtric"){
  
  enrichment_results<-list()
  
  all_genes <-unique(unlist(pathways))
  
  for (pathway_name in names(pathways)){
    
    list_of_genes <- pathways[[pathway_name]]
    
    overlap <- intersect(genes_of_interest, list_of_genes)
    
    #calcualte enrichment statistic usign hypergeometric test
    
    if (test_type == "hypergeometric"){
      
      p_value <- phyper(length(overlap) - 1, length(all_genes), 
                        length(list_of_genes), length(genes_of_interest), lower.tail = FALSE)
      
    }else{
      
      fisher_result <- fisher.test(matrix(c(length(overlap), 
                                            length(setdiff(genes_of_interest, overlap)),
                                            length(setdiff(list_of_genes, overlap)),
                                            length(setdiff(union(genes_of_interest, 
                                                                 list_of_genes), overlap))),
                                          nrow = 2))
      
      p_value = fisher_result$p.value
      
    }
    
    
    # Store results
    enrichment_results[[pathway_name]] <- list(
      list_of_genes = list_of_genes,
      overlap_genes = overlap,
      p_value = p_value
    )
    
  }
  
  return(enrichment_results)
}

clean_enrichment<-function(list_of_genes, enrichment_results, list_of_interest, grouped = FALSE, correction  = "bonferroni"){
  
  p_values<-sapply(enrichment_results, function(result) result$p_value)
  
  adjusted_p_values <- p.adjust(p_values, method = correction)  # Bonferroni correction
  
  adjusted_p_values<-as.data.frame(adjusted_p_values)
  colnames(adjusted_p_values)[1]<-"p_vals"
  adjusted_p_values$pathways <- rownames(adjusted_p_values)
  
  adjusted_p_values<-adjusted_p_values[order(adjusted_p_values$p_vals),]
  
  top_pathway = list()
  
  #filtered_list<-unlist(filtered_list, recursive = FALSE)
  
  chosen_list<-list_of_interest
  
  if (grouped == FALSE){
    for (i in 1:length(list_of_genes)){
      
      for (j in 1:nrow(adjusted_p_values)){
        
        tmp_path = rownames(adjusted_p_values)[[j]]
        
        if (as.character(list_of_genes[[i]]) %in% chosen_list[[tmp_path]]){
          
          top_pathway[[i]] = unlist(list(list_of_genes[[i]], adjusted_p_values[j, "pathways"][[1]]))
          
          break
          
        }
        
      }
    }
    
    top_pathway_df<-do.call(rbind, top_pathway)
    
    return(top_pathway_df)
    
  }else{
    
    # Find names of pathways that contain all genes from your group
    
    for (j in 1:nrow(adjusted_p_values)){
      
      tmp_path = rownames(adjusted_p_values)[[j]]
      
      if (all(list_of_genes %in% list_of_interest[[tmp_path]])){
        print("yes")
        
        for (k in 1:length(list_of_genes)){
          top_pathway[[k]]<-unlist(list(list_of_genes[[k]], tmp_path))
        }
        
        break
        
      }else{
      }
      
    }
    # View the result
    
    return(top_pathway)
    
    
    
  }
  
}

geneid_to_protein <- function(gene_ids, mapping_list) {
  # Initialize an empty vector to store protein names
  protein_names <- character(length(gene_ids))
  
  # Iterate over each gene ID and find the corresponding protein name
  for (i in seq_along(gene_ids)) {
    # Use match to find the index of the gene ID in the mapping list
    match_index <- match(gene_ids[i], unlist(mapping_list))
    
    # Check if a match was found
    if (!is.na(match_index)) {
      # Use the match index to get the corresponding protein name
      protein_names[i] <- names(mapping_list)[match_index]
    } else {
      # If no match was found, set the protein name to NA
      protein_names[i] <- NA
    }
  }
  
  return(protein_names)
}

name_check<-function(word){
  tmp_word<-gsub("-","", word)
  tmp_word<-gsub("_","", tmp_word)
  
  replacement<-updated_names[updated_names$Current==tmp_word, "Update"]
  
  return(replacement)
}


