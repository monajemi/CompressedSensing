.packageName <- 'PT'

predictPT <-function(delta, field, type="w"){

#PredPT <- data(PredPT);
if(type=="w"){

if(field=="R"){
x = PredPT$delta.Real
y = PredPT$eps.Real
# fit a spline to data of (delta,epsilon)
f<-splinefun(x=x,y=y) 
# predict at new delta
}


if(field=="C"){
x = PredPT$delta.Cplex
y = PredPT$eps.Cplex
# fit a spline to data of (delta,epsilon)
f<-splinefun(x=x,y=y) 
# predict at new delta
}

if(field=="Q" || field=="H2"){
x = PredPT$delta.Q
y = PredPT$eps.Q
# fit a spline to data of (delta,epsilon)
f<-splinefun(x=x,y=y) 

}

if(field=="O" || field=="H3"){
    x = PredPT$delta.O
    y = PredPT$eps.O
    # fit a spline to data of (delta,epsilon)
    f<-splinefun(x=x,y=y)
    # predict at new delta
}

if(field=="Bnd"){
f<- function(x){2*x-1}

}

if(field=="Pos"){
x = PredPT$delta.Pos
y = PredPT$eps.Pos
# fit a spline to data of (delta,epsilon)
f<-splinefun(x=x,y=y) 

}


}else if(type=="s"){

if(field=="R"){
    x = PredPT$delta.Real.S
    y = PredPT$eps.Real.S
    # fit a spline to data of (delta,epsilon)
    f<-splinefun(x=x,y=y)
    # predict at new delta
}


}

u = f(delta);

return(u * (u > 0) )

}
