
###################################################
### Combine STRING scores
###################################################
#script escrito por Iara Souza
# iaradsouza1@gmail.com 

combinescores<-function(dat, evidences="all", confLevel=0.4){
      if(evidences[1]=="all"){
            edat<-dat[,-c(1,2,ncol(dat))]
      } else {
            if(!all(evidences%in%colnames(dat))){
                  stop("NOTE: one or more 'evidences' not listed in 'dat' colnames!")
            }
            edat<-dat[,evidences]
      }
      edat<-edat/1000
      edat<-1-edat
      sc<- apply(X = edat, MARGIN = 1, FUN = function(x) 1-prod(x))
      dat<-cbind(dat[,c(1,2)],combined_score=sc)
      idx<-dat$combined_score>=confLevel
      dat<-dat[idx,]
      return(dat)
}




