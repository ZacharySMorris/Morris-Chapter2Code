###Centered plot for slope comparisons###

MA_number <- 1
ConfIntList <- Clades3.MajorAxis
LineColor <- Combined_Morphospace_Clades3$legendcolor

for (i in 1:length(ConfIntList$Groups)){

plot(0, 0, type = "n",
     xlim = c(-0.3,0.3),
     ylim = c(-0.3,0.3),
     xlab = "",
     ylab = "",
     axes = FALSE,
     frame.plot = FALSE,
     asp=F)

abline(h=0, lty=3)
abline(v=0, lty=3)

Adult_mshape <- mshape(Adult_CrocDorsal.pca$pc.scores)

polygon(c(rev(Adult_Total_Morphospace.MajorAxis_1$NEW_Xs$Crocodylia-Adult_mshape[1]),
          Adult_Total_Morphospace.MajorAxis_1$NEW_Xs$Crocodylia-Adult_mshape[1]),
        c(rev(Adult_Total_Morphospace.MajorAxis_1$preds_CIupper$Crocodylia-Adult_mshape[2]),
          Adult_Total_Morphospace.MajorAxis_1$preds_CIlower$Crocodylia-Adult_mshape[2]),
        col = alpha('grey',0.4), border = NA)

lines(Adult_Total_Morphospace.MajorAxis_1$NEW_Xs$Crocodylia-Adult_mshape[1],
      (Adult_Total_Morphospace.MajorAxis_1$preds_CIupper$Crocodylia-Adult_mshape[2]),
      lty = 'dashed',
      col = unique(Maturity_combined$color$Adult))
lines(Adult_Total_Morphospace.MajorAxis_1$NEW_Xs$Crocodylia-Adult_mshape[1],
      (Adult_Total_Morphospace.MajorAxis_1$preds_CIlower$Crocodylia-Adult_mshape[2]),
      lty = 'dashed',
      col = unique(Maturity_combined$color$Adult))
abline(0,Adult_Total_Morphospace.MajorAxis_1$Slope_list$Crocodylia[MA_number],
       col = unique(Maturity_combined$color$Adult), lwd=3, lty=1)

Ontogenetic_Ecomorph_list <- c("Blunt","Moderate","Slender","Tomistoma","Gavialis")

for(j in 1:length(Ontogenetic_Ecomorph_list)){
  current_group <- Ontogenetic_Ecomorph_list[j]
  current_mshape <- mshape(Crocodylia_Ecomorphs_G_T_Unique$PCvalues[[current_group]])


  polygon(c(rev(Crocodylia_Ecomorphs_G_T_Unique_Total_Morphospace.MajorAxis_1$NEW_Xs[[current_group]]-current_mshape[1]),
            Crocodylia_Ecomorphs_G_T_Unique_Total_Morphospace.MajorAxis_1$NEW_Xs[[current_group]]-current_mshape[1]),
          c(rev(Crocodylia_Ecomorphs_G_T_Unique_Total_Morphospace.MajorAxis_1$preds_CIupper[[current_group]]-current_mshape[2]),
            Crocodylia_Ecomorphs_G_T_Unique_Total_Morphospace.MajorAxis_1$preds_CIlower[[current_group]]-current_mshape[2]),
          col = alpha('grey',0.4), border = NA)

  lines(Crocodylia_Ecomorphs_G_T_Unique_Total_Morphospace.MajorAxis_1$NEW_Xs[[current_group]]-current_mshape[1],
        (Crocodylia_Ecomorphs_G_T_Unique_Total_Morphospace.MajorAxis_1$preds_CIupper[[current_group]]-current_mshape[2]),
        lty = 'dashed',
        col = unique(Crocodylia_Ecomorphs_G_T_Unique$color[[current_group]]))
  lines(Crocodylia_Ecomorphs_G_T_Unique_Total_Morphospace.MajorAxis_1$NEW_Xs[[current_group]]-current_mshape[1],
        (Crocodylia_Ecomorphs_G_T_Unique_Total_Morphospace.MajorAxis_1$preds_CIlower[[current_group]]-current_mshape[2]),
        lty = 'dashed',
        col = unique(Crocodylia_Ecomorphs_G_T_Unique$color[[current_group]]))
  abline(0,Crocodylia_Ecomorphs_G_T_Unique_Total_Morphospace.MajorAxis_1$Slope_list[[current_group]][MA_number],
         col = unique(Crocodylia_Ecomorphs_G_T_Unique$color[[current_group]]), lwd=3, lty=1)
}

  lengh_x <- max(ConfIntList$NEW_Xs[[i]]) - min(ConfIntList$NEW_Xs[[i]])
  clip(-(lengh_x/2),lengh_x/2,-0.3,0.3)

  abline(0,ConfIntList$Slope_list[[i]][MA_number], col = LineColor[i], lwd=5, lty=1 )
  mtext(ConfIntList$Groups[[i]], at = -0.25)
}



plot(0, 0, type = "n",
     xlim = c(-0.3,0.3),
     ylim = c(-0.3,0.3),
     xlab = "",
     ylab = "",
     axes = FALSE,
     frame.plot = FALSE,
     asp=F)

abline(h=0, lty=3)
abline(v=0, lty=3)
