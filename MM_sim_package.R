
# for n = 4
for(i in 6:10){
M0 = c(0,0.3,0.5,1,5)
  
  SimMdata(n = 4, parsN = c(1,1,0,0), age = 30, M0= M0[i-5],num = i)

}


# for n = 3
for(i in 6:10){
  M0 = c(0,0.3,0.5,1,5)
  
  SimMdata(n = 3, parsN = c(1,1,0), age = 30, M0= M0[i-5],num = i)
  
}


# for n = 2
for(i in 5:1){
  M0 = c(0,0.3,0.5,1,5)
  
  SimMdata(n = 2, parsN = c(1,1), age = 30, M0= M0[6-i],num = i)
  
}
