

list_of_matrix<-list(matrix(1:1000,nrow=100),matrix(1:1000,nrow=100))
for (matrix in list_of_matrix){
  print (matrix)
  #do some other task
  print(matrix*2000)}
