
library(performanceEstimation)
suppressWarnings(suppressMessages(library(uba)))

#nombre de los ficheros
f = "y.txt"
f_pred = "y_pred.txt"
out = "phi.txt"
f1_out = "f1.txt"
th.rel = 0.8


#Leemos las etiquetas y las predichas:
y = read.table( file = f, sep = "")[,1]
y_pred = read.table( file = f_pred, sep = "")[,1]

#Calculamos parametros de la funci?n de utilidad, relevancia y de p?rdida:
pc = phi.control(y, method="extremes")
lc <- loss.control(y) 

#Calculamos la funci?n de relevancia
#phi_y = phi(y, pc)

#Calculamos el F1 score:
#Calculamos la funci?n de utilidad:
up <- util.control(umetric="Fm",event.thr=th.rel,beta=1,p=0.5) 
f1_score = util(y_pred,y,phi.parms=pc,loss.parms=lc,util.parms = up) 





# if (file.exists(f)){
#   #Delete file if it exists
#   file.remove(f)
#   file.remove(f_pred)
# }
borrar = function(file1, file2){
  file.remove(file1)
  file.remove(file2)
  
  1
}

trash=ifelse(file.exists(f),borrar(f, f_pred), 1)

#Write phi output:
#write.table(phi_y,row.names=FALSE, col.names=FALSE, file = out, sep = "")
write.table(as.data.frame(f1_score),row.names=FALSE, col.names=FALSE, file = f1_out, sep = "")

