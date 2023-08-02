## FGM Copula
## Case1 + Cox + Copula

library(survival) 
library(MASS)  
library(pracma)  
library(nleqslv)
library(statmod) 
library(splines2) 

alln = c(50, 100, 200, 400, 10^6, 10, 20)  
allbeta0 = c(0, 0.4) 
allgamma0 = c(0, 0.4, -0.4) 


alltau0 = c(0.1, 0.2, -0.1, -0.2, 0)  ## 


#####################################

n = alln[4]


mj = 1000  




#################

for(iii in 1:5)
{
  for(II in 1:2)
  {
    JJ = 1
    
      
      tau0 = alltau0[iii]
      alpha0 = 4.5*tau0

      

      beta0 = allbeta0[II]
      

      gamma0 = allgamma0[JJ]
      
      
      
      
      es = matrix(0,mj,2)
      ev = matrix(0,mj,2)
      cover = matrix(0,mj,2)
      
      
      # fgm copula
      copula<-function(u,v,af){
        u*v + af*u*v*(1-u)*(1-v)
      }
      

      mal<-function(u,v,af){
        u + af*u*(1-u)*(1-2*v) 
      }
      
      
      
      for(ct in 1:mj)
      {
        
        Z = numeric(n)
        
        TT = numeric(n)  
        A = numeric(n)  
        C = numeric(n)  
        
        
       
        
        nex = 2*n  
        
        lok = 0  
        
        
        kk=0.35
        
        bb=0.05
        
        
        
        ll1 = 0
        
        while(lok<n)
        {
          Z0ex = rbinom(nex, 1, 0.5)   
          
          FT = runif(nex)
          T0ex = -log(1-FT)/exp(beta0*Z0ex)/kk
          
          UC = runif(nex)
          Au = 1 + alpha0*(1-2*FT)
          Bu = sqrt( Au^2 - 4*UC*(Au-1) )
          FC = 2*UC/(Au+Bu)
          C0ex = (-log(1-FC)/exp(gamma0*Z0ex)/bb)^(1/2)
 
          
          A0ex = rexp(nex, rate=1)
          
          TBA = (T0ex>=A0ex)*(C0ex>A0ex)
          LTBA = sum(TBA)
          
          if(LTBA+lok>=n){
            TT[(lok+1):n] = T0ex[TBA==1][1:(n-lok)]
            A[(lok+1):n] = A0ex[TBA==1][1:(n-lok)]
            C[(lok+1):n] = C0ex[TBA==1][1:(n-lok)]
            Z[(lok+1):n] = Z0ex[TBA==1][1:(n-lok)]
            lok = n
          }else{
            TT[(lok+1):(lok+LTBA)] = T0ex[TBA==1]
            A[(lok+1):(lok+LTBA)] = A0ex[TBA==1]
            C[(lok+1):(lok+LTBA)] = C0ex[TBA==1]
            Z[(lok+1):(lok+LTBA)] = Z0ex[TBA==1]
            lok = lok+LTBA
          }
          
          ll1=ll1+1

        }
        
        
        

        delta = as.numeric(TT<=C)
        
        

        data.AC = sort(unique(c(A,C)))
        
        

        mknots = 3
        

        

        knots.all0 = seq( from=min(data.AC), to=max(data.AC), by=(max(data.AC)-min(data.AC))/(mknots+1) )

        knots.all1 = knots.all0[c(-1,-length(knots.all0))]
        
        

        knots.all = quantile(data.AC, probs=seq(from=1/(mknots+1), by=1/(mknots+1), length=mknots) )
        
        

        degree.em = 2
        
        

        Sbasis.all = iSpline(c(A,C), knots = knots.all, degree = degree.em, intercept = TRUE)


        dv.Sbasis.all = iSpline(c(A,C), knots = knots.all, degree = degree.em, intercept=TRUE, derivs=1)

        
        
        

        Sbasis.C = Sbasis.all[(n+1):(2*n),]
        
        dv.Sbasis.C = dv.Sbasis.all[(n+1):(2*n),]
        

        Sbasis.A = Sbasis.all[1:n,]
        
        
        co = function(v)   
        {
          v^2
        }
        
        
        ##
        LT0.C = function(v)
        {
          as.numeric( Sbasis.C%*%co(v) )
        }
        
        ##
        LT0.A = function(v)
        {
          as.numeric( Sbasis.A%*%co(v) )
        }
        
        
        ##
        LC0.C = function(u)
        {
          as.numeric( Sbasis.C%*%co(u) )
        }
        
        ##
        LC0.A = function(u)
        {
          as.numeric( Sbasis.A%*%co(u) )
        }
        
        
        

        dlc0.C = function(u)
        {
          as.numeric( dv.Sbasis.C%*%co(u) )
        }
        
        

        lcosp = degree.em + mknots + 1

        
        
        
        neg.ll = function(pat)
        {
          ade = 1e-7
          
          th = pat[1:2]
          v = pat[3:(lcosp+2)]
          u = pat[(lcosp+3):(2*lcosp+2)]
          
          bt = th[1]
          gm = th[2]
          
          exp.bt.Z = exp(bt*Z)
          exp.gm.Z = exp(gm*Z)
          
          FT.C = 1 - exp( -LT0.C(v)*exp.bt.Z )
          FT.C = pmax( ade, pmin(FT.C, 1-ade) )
          
          FT.A = 1 - exp( -LT0.A(v)*exp.bt.Z )
          FT.A = pmax( ade, pmin(FT.A, 1-ade) )
          
          FC.C = 1 - exp( -LC0.C(u)*exp.gm.Z )
          FC.C = pmax( ade, pmin(FC.C, 1-ade) )
          
          FC.A = 1 - exp( -LC0.A(u)*exp.gm.Z )
          FC.A = pmax( ade, pmin(FC.A, 1-ade) )
          

          fc = dlc0.C(u)*exp.gm.Z*(1-FC.C)
          fc = pmax(0, fc)
          
          mcc = mal(FT.C, FC.C, alpha0)
          mcc = pmax( 0, pmin(mcc,1) )
          
          mac = mal(FT.A, FC.C, alpha0)
          mac = pmax( 0, pmin(mac,1) )
          
          caa = copula(FT.A, FC.A, alpha0)
          caa = pmax( 0, pmin(caa,1) )
          
          
          ll = sum(  log(fc) + delta*log(mcc-mac+1e-320) 
                     + (1-delta)*log(1-mcc) - log(1-FT.A-FC.A+caa)  )
          
          -ll
        }
        
        
        neglgm = function(pg)
        {
          gm = pg[1]
          u = pg[2:(lcosp+1)]
          
          gm.Z = gm*Z
          
          lgm = sum( log(dlc0.C(u)) + gm.Z + exp(gm.Z)*(LC0.A(u)-LC0.C(u)) )
          
          -lgm
        }
        
        
        
        neglbt = function(pb)
        {
          ade = 1e-7
          
          bt = pb[1]
          v = pb[2:(lcosp+1)]
          
          exp.bt.Z = exp(bt*Z)
          LT0Av = LT0.A(v)
          LT0Cv = LT0.C(v)
          
          ST.A = exp( -LT0Av*exp.bt.Z )
          ST.A = pmax( 0, pmin(ST.A, 1) )
          
          ST.C = exp( -LT0Cv*exp.bt.Z )
          ST.C = pmax( 0, pmin(ST.C, 1) )
          
          lbt = sum(  delta*( log(ST.A-ST.C+1e-120) + LT0Av*exp.bt.Z )  
                      + (1-delta)*exp.bt.Z*( LT0Av-LT0Cv )  )
          
          -lbt
        }
  

        min.u = function(u)
        {
          sum( (bb*C^2 - as.numeric(Sbasis.C%*%co(u)))^2 )
        }
        
        

        u00 = optim(par=rep(0.1, lcosp), min.u, method="BFGS")$par
        

        my.surv.C = Surv(A, C, rep(1,n), type='counting')  

        coxph.fit.C = coxph( my.surv.C~Z )
        gm00 = coxph.fit.C$coef  
        
        gm.u00 = c(gm00, u00)
        ini.gm.u = nlm(neglgm, gm.u00)$estimate
        
        ##################################################################

        ini.gm = ini.gm.u[1]
        ini.u = ini.gm.u[2:(lcosp+1)]
        ##################################################################
        
        
        
        

        min.v = function(v)
        {
          sum( (kk*C - as.numeric(Sbasis.C%*%co(v)))^2 )
        }
        
        

        v00 = optim(par=rep(0.1, lcosp), min.v, method="BFGS")$par
        
        bt.v00 = c(beta0, v00)

        ini.bt.v = optim(bt.v00, neglbt)$par
        
        ##################################################################

        ini.bt = ini.bt.v[1]
        ini.v = ini.bt.v[2:(lcosp+1)]
        ##################################################################
        

        
        allpar0 = c(ini.bt, ini.gm, ini.v, ini.u)

        
        es.crude.all = optim(allpar0, neg.ll)$par
        es.crude = es.crude.all[1:2]
        btcr = es.crude[1]
        gmcr = es.crude[2]
        
        if(  (btcr<(beta0-5)) | (btcr>(beta0+5)) 
             | (gmcr<(gamma0-5)) | (gmcr>(gamma0+5))  ){
          es[ct,] = c(500, 500)
          ev[ct,] = c(NaN, NaN)
        }else{
          es[ct,] = es.crude
          ifm = hessian(neg.ll, es.crude.all)
          inver.ifm = ginv(ifm)[1:2, 1:2]
          ev[ct,] = sqrt(diag(inver.ifm))
        }
        
        

      }
      
      
      
      ores = c(beta0, gamma0)
      oresma = rep(1,mj)%*%t(ores)

      
      
      esmean = as.numeric(rep(1,mj)%*%es)/mj
      
      bias = esmean - ores
      SD = sqrt( as.numeric((rep(1,mj)%*%(es^2))-mj*esmean^2)/mj )
      ASE = as.numeric(rep(1,mj)%*%ev)/mj
      
      qu = qnorm(0.975)
      leftma = es - qu*ev
      rightma = es + qu*ev
      cover = (leftma<=oresma)*(oresma<=rightma)
      
      CP = as.numeric(rep(1,mj)%*%cover)/mj
      
   
      
      error1 = unique(which(is.na(ev)==TRUE, arr.ind=TRUE)[,1])

      
      error.bt = which( (es[,1]<(beta0-1.5)) | (es[,1]>(beta0+1.5))  )
      error.gm = which( (es[,2]<(gamma0-1.5)) | (es[,2]>(gamma0+1.5)) )
      
      
      
      error = unique( sort(c(error1, error.bt, error.gm)) )
      
      
      
      if( length(error)!=0 ) { 
        bias = apply(es[-error,], 2, mean) - ores
        SD = apply(es[-error,], 2, sd)
        ASE = apply(ev[-error,], 2, mean)
        CP = apply(cover[-error,], 2, mean) }
      
      
      
      bias = round( bias, digits=4 )
      SD = round( SD, digits=4 )
      ASE = round( ASE, digits=4 )
      CP = round( CP, digits=4 )
      
      
      
      name=paste('Cox-FGM', 'n=', n, 'bata=', beta0, 'gamma=', gamma0, 'tau=', tau0, '.txt')
      re00 = cbind(bias, SD, ASE, CP)
      result = t( c(re00[1,], re00[2,]) )
      write.table(result, name)
      
      
    }
    
  }
  

