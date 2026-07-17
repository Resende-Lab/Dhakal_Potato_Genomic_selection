
#########################<-Potato breeding program
P = runif(6) #p-values for GxY effect

for(year in 1:6){

cat("FillPipeline year:",year,"of 6\n")
  ####<----------------YEAR 1: Cross
  Clones_C0 <- randCross(pop = Parents, nCrosses= 1000, nProgeny = 20, ignoreSexes = TRUE)
  # Mass selection
  Clones <- selectInd(Clones_C0, 13000, use = "rand")
  

if(year<6){
  p = P[6-year]
  ####<---------------YEAR 2: Production Seedlings and First Clonal
  Field_seedling <- setPheno(Clones, h2 = h2, reps = 5)
  Field_seedling_Sel <- selectInd(Field_seedling, 1000, use = "pheno",
                                  trait = selIndex, b = weight, scale = TRUE)
}

if(year<5){
  p = P[5-year]
  ####<---------------YEAR 3: Second clonal generation
  First_Clonal <- setPheno(Field_seedling_Sel, reps = 2)
  First_Clonal_Sel <- selectInd(First_Clonal,  400, use = "pheno",
                                trait = selIndex, b = weight, scale = TRUE)
  
}

if(year<4){
  p = P[4-year]
  ####<---------------YEAR 4: Third clonal generation
  Second_Clonal <- setPheno(First_Clonal_Sel, reps = 6)
  Second_Clonal_Sel <- selectInd(Second_Clonal, 90, use = "pheno",
                                 trait = selIndex, b = weight, scale = TRUE)
}

if(year<3){
  p = P[3-year]
  ####<---------------YEAR 5: Official trials 1
  Third_Clonal <- setPheno(Second_Clonal_Sel, reps = 15)
  Third_Clonal_Sel <- selectInd(Third_Clonal, 15, use = "pheno",
                                trait = selIndex, b = weight, scale = TRUE)
}


if(year<2){
  p = P[2-year]
  ####<---------------YEAR 6: Official trials 2
  Fouth_Clonal <- setPheno(Third_Clonal_Sel, reps = 30)
  Fouth_Clonal_Sel <- selectInd(Fouth_Clonal, 5, use = "pheno",
                                trait = selIndex, b = weight, scale = TRUE)
  
}

#Year 7, release

}
