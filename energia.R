#### OJO CON LAS CABECERAS RARAS #####
#puede ser que por una cabecera extraña, como el asignar automaticamente
#el nombre raro de un vector a una columna genere que no pueda leer el archivo.
#Ojo con los simbolos y los espacios.

#Script para contabilizar amino acidos



exe01= ("sed 's/Dapma.*//' dmagset7finloc9b.puban.aa > aa1.fasta")
exe02=("cat aa1.fasta | sed ':a;N;$!ba;s/\n//g'> aa2.fasta")
exe03=("sed 's/>/\n/g' aa2.fasta > aa3.fasta")

#Conteo por línea
exe04=("awk -F'|' 'BEGIN{print 'count' 'lineNum'}{print gsub(/A/,'' '\t' NR}' aa3.fasta > peo.csv")
system(exe04)
#Eliminas la primera y segunda linea (linea cabecera y  linea espacio).

exe05=("sed -i '1,2d' Aa_A.csv")

#Hago un cat de todos los archivos generados

#No se puede hacer este cat !!! porque no se sabe como ordena.
exe06=("cat Aa* > cat_aa.csv")

#Debe ser uno por uno en el orden que uno defina
#Se hizo en esta carpeta para hacer el cat mas limpio ~/Escritorio/Mauricio/PB/cat

exe07=("cat Aa_A.csv Aa_C.csv > AC.csv")


#cat final que incluye a todos:    conteo.csv




# Me quedo solo con la primera columna que es el conteo.

exe08=("awk '{print$1}' conteo.csv > conteo1.csv")



#******


#CONTERO DE ATP POR SECUENCIA
#Es un conteo de genes totales


setwd("~/Desktop/Mauricio/PB_SE/PB_SE")
#conteo<-read.table("salida_conteo_todos_aa.txt", header = F)
conteo<-read.csv("conteo1.csv", header = F)
#conteo2=as.character(conteo1)
#conteo=data.frame(conteo2)

#El vector conteo contiene el conteo de cada uno de los aminácidos para todos los
#genes, no separado por sobre y sub.
#Esto fue obtenido con el srcipt de perl "contador_aa.pl" de la misma carpeta


#View(conteo)

A=c(conteo[1:29127,1])
#View(C)

C=c(conteo[29128:58254,1])

R=c(conteo[58255:87381,1])

N=c(conteo[87382:116508,1])

D=c(conteo[116509:145635,1])

Q=c(conteo[145636:174762,1])

E=c(conteo[174763:203889,1])

G=c(conteo[203890:233016,1])

H=c(conteo[233017:262143,1])

I=c(conteo[262144:291270,1])

L=c(conteo[291271:320397,1])

K=c(conteo[320398:349524,1])

M=c(conteo[349525:378651,1])

f=c(conteo[378652:407778,1])

P=c(conteo[407779:436905,1])

S=c(conteo[436906:466032,1])

t=c(conteo[466033:495159,1])

W=c(conteo[495160:524286,1])

Y=c(conteo[524287:553413,1])

V=c(conteo[553414:582540,1])

##Calculo gasto por AA
#df=(cbind(C,A,R,N,D,Q,E,G,H,I,L,K,M,f,P,S,t,W,Y,V))
#Los números corresponden al número de moléculas necesarias para sintetizar el Aa

df2=(data.frame(C*24,A*12.5,N*4,D*1,Q*8.5,E*9.5,G*14.5,P*12.5,S*15,Y*56.5))



#Dos citas para aminoacidos escenciales
#doi.org/10.1007/BF02026769
#www.fao.org/3/AB492S/AB492S01.htm
#Elimino:
#Treonina Valina Leucina Isoleucina Metionina Triptófano
#Lisina Histidina Arginina Fenilalaina




#ACA OBTENTO LA TABLA CON LOS VALORES TOTALES DE costo de todos los AA POR GEN,
#SIN LA EXPRESION.
sss<-rowSums(df2)

#View(sss)

conteoscd1<-read.table("pb.csv", header = T, sep = ",", quote = NULL)
#View(conteoscd1)
#Lo convierto a data frame ya que sss no tiene header, por lo que no hace el cbind
#ya que no estarían valanceados. Por eso el dataframe con el header =T para sss.
sss1=as.data.frame(sss, header = F)


p=cbind(conteoscd1[,1],sss1)
#View(p)

write.csv(p, file="unitario_pb_energia_de.csv", quote = T, row.names = F)

#Ahora busco los datos de expresión normalizada



library("DESeq2")
# Leemos conteos
conteoscd<-read.csv("pb.csv", header = T, row.names = "GEN", sep = ",")


#View(conteoscd)
dim(conteoscd)# dimensiones de la tabla

# La serie "X1000_follow_up" son los No Meditadores, y
# la serie "X2000_follow_up" son los Meditadores regulares.
colnames(conteoscd[1])  #nombre de columna 1
colnames(conteoscd[2]) #nombre de columna 56
colnames(conteoscd[7]) #nombre de columna 57
#col 2 a 56 son los no meditadores 1 mil ...
#57 adelante meditadores 2 mil..


#Creaci?n de objeto countdata, a partir del Subseting de los datos "ConteosMeditacion".
# Se selecciona Todas las filas y desde la columna 2 hasta la 7:
countdata<-conteoscd[, c(1:6)]
#filtrado de conteos
#me quedo solo con los conteos, no con los nombres
dim(countdata) #dimensiones del objeto "countdata"


#coldata: tabla con informacion de las muestras
coldata = data.frame(
  #crear data frame con dos columnas
  row.names = colnames( countdata ),
  #nombres con lo anterior
  tipo = c( rep("Zontrol", 3), rep("levels", 3)))
# a los datos de countdata, escribir No meditador 55 veces, etc
# View(coldata)

# Construimos el objeto DESeqDataSet a partir de la matriz de conteos y la informaci?n de las muestras
dds<- DESeqDataSetFromMatrix(countData = countdata,
                             #Entrega conteos, metadata, y permite comparar entre tipo Nomed y medit (que es desing=tipo)
                             colData = coldata,
                             design = ~ tipo)
dds=DESeqDataSetFromMatrix(counts,targets,design=~tipo)
t <- proc.time() # Inicia el cron?metro
dds <- DESeq(dds)
proc.time()-t    # Detiene el cron?metro.
#En mi compu demora alrededor de 74.42 s aprox.



##########
#OBTENCION DE DATOS NORMALIZADOS - OPCIONAL
##########



# Obtener conteos normalizados (opcional, si se necesita el dato de normalizacion)
sf <- estimateSizeFactors(dds)
normctn<-counts(sf, normalized=TRUE)
#View(normctn)
#Si imprimo aca conserva los nombres.
write.table(normctn, file="tratamiento_normalizados_PB.txt", quote = T, sep="\t", row.names = T)


#####MULTIPLICACION
#LEO EL ARCHIVO ANTERIOR

norm_pb=read.table("tratamiento_normalizados_PB.txt", header = T)
#View(norm_pb)

#multiplicamos
#sss = Es el valor del costo de ATP total por gen
#norm_pb = Valor de la expresion normalizada
energia_de= norm_pb*sss
#View(energia_de)
#energia_de = costo de expresion con nombres
#                    XI_CO_r1_04 XI_CO_r2_05 XI_CO_r3_06 XI_PB_r1_07 XI_PB_r2_08 XI_PB_r3_09
#Dapma7bEVm000001t1 10998910.8 12984174.0 13821879.8 11785035.12 10540504.77 10102120.88

#Para tener solo los nombres de que genes están solo en control y solo en Pb debo:
#Calcular logaritmo
#Restarle -1
#Calcular logaritmo
#eliminar todas las filas que tengan NA y INF (0 no)

log_energia_de=log(energia_de)
uno_log_energia_de=log_energia_de-1
na_uno_log_energia_de=log(uno_log_energia_de)
#View(na_uno_log_energia_de)
#Ahora los separo en dos grupos para que al eliminar -Inf y Na por separado me queden
#Control y tratamiento con sus correspondientes genes expr. que pueden no estar en
#El otro.
control=as.data.frame(na_uno_log_energia_de[,1:3])
#View(control) #costo ponderado de expresion de solo control
write.table(control, file="control_na_uno_log_energia_de.txt", quote = T, sep="\t", row.names = T)
#Ahora en linux elimino las filas que contengan NA y -Inf
exe = ("sed '/-Inf/d; /NA/d' control_na_uno_log_energia_de.txt > control_sin_na_uno_log_energia_de.txt")
system(exe)
#control_sin_na_uno_log_energia_de.txt = cotiene todas las celdas menos las que tienen -Inf y NA
#Nombre y costo de 3 replicas

pb=as.data.frame(na_uno_log_energia_de[,4:6])
#View(pb) #costo ponderado de expresion de solo pb
write.table(pb, file="pb_na_uno_log_energia_de.txt", quote = T, sep="\t", row.names = T)
#Ahora en linux elimino las filas que contengan NA y -Inf
exe1=("sed '/-Inf/d; /NA/d' pb_na_uno_log_energia_de.txt > pb_sin_na_uno_log_energia_de.txt")
system(exe1)
#pb_sin_na_uno_log_energia_de.txt = cotiene todas las celdas menos las que tienen -Inf y NA
#Nombre y costo de 3 replicas


#############             VENN           ####################

#Hago un grep entre solo los nombres de los archivos obtenidos y su valor normalizado.
#Con esto obtendo los nombres de los genes expresados en cada condión
exe2=("awk '{print $1}' control_sin_na_uno_log_energia_de.txt > nombres_solo_control")
system(exe2)
exe3=("awk '{print $1}' pb_sin_na_uno_log_energia_de.txt > nombres_solo_pb")
system(exe3)

exe6 =("grep -wf nombres_solo_control tratamiento_normalizados_PB.txt > solo_control_norm.csv")
system(exe6)
exe7= ("grep -wf nombres_solo_pb tratamiento_normalizados_PB.txt > solo_pb_norm.csv")
system(exe7)





solo_control_norm=read.csv("solo_control_norm.csv", header = T,sep = "\t")
#View(solo_control_norm)
solo_pb_norm=read.table("solo_pb_norm.csv", header = T,sep = "\t")
#View(solo_pb_norm)




###LOS ARCHIVOS ANTERIORES TIENEN QUE MULTIPLICARSE POR EL VALOR UNITARIO DE ENERGIA

#Obtener los valores unitarios de expresion por gen para cada tratamiento
#Luego de esto hay que multiplicar
#hacemos un grep con los nombres de cada trat contra los valores de energia uniq

exe4=("grep -wf nombres_solo_control unitario_pb_energia_de.csv > control_de_en.csv")
system(exe4)
exe5=("grep -wf nombres_solo_pb unitario_pb_energia_de.csv > pb_de_en.csv")
system(exe5)

#Leo los archivos anteriores del grep

control_en_unit=read.csv("control_de_en.csv", header = F, dec =".")
control_en_unit_1=as.data.frame(control_en_unit)
control_en_unit_1[, 2]  <- as.numeric(control_en_unit_1[, 2])


#Costo unitario de expresion de genes en Co
#View(control_en_unit)
pb_en_unit=read.csv("pb_de_en.csv", header = F, sep = ",")

#Costo unitario de expresion de genes en Pb
#View(pb_en_unit)


########      Costo por la EXPRESION DIFERENCIAL     #######

control1=control_en_unit[,2]*solo_control_norm[,1:3]
control1_1=control_en_unit[,2]
#View(control1)

pb_1=pb_en_unit[,2]*solo_pb_norm[,4:6]
pb_1_1=pb_en_unit[,2]
#View(pb_1_1)




#vector_de_pb=c(pb_1$XI_PB_r1_07,pb_1$XI_PB_r2_08,pb_1$XI_PB_r3_09)
vector_de_pb_1=c(pb_1_1)
#View(vector_de_pb)
vector_de_pb=c(pb_1[,1],pb_1[,2],pb_1[,3])
#View(vector_de_pbb)

vector_de_pb2=as.data.frame(vector_de_pb)
vector_de_pb2_1=as.data.frame(vector_de_pb_1)
#View(vector_de_co)
conteo_vector_de_pb=nrow(vector_de_pb2)
conteo_vector_de_pb_1=nrow(vector_de_pb2_1)
conteo_vector_de_pb
conteo_vector_de_pb_1
vector_de_pb1=round(vector_de_pb, digits = 3)
vector_de_pb1_1=round(vector_de_pb_1, digits = 3)


vector_de_co=c(control1[,1],control1[,2],control1[,3])
vector_de_co_1=c(control1_1)

#Lo convierto a data.frame sino no lo va a poder contar.
vector_de_co2=as.data.frame(vector_de_co)
vector_de_co2_1=as.data.frame(vector_de_co_1)
#View(vector_de_co)
conteo_vector_de_co=nrow(vector_de_co2)
conteo_vector_de_co_1=nrow(vector_de_co2_1)
conteo_vector_de_co
conteo_vector_de_co_1
vector_de_co1=round(vector_de_co, digits = 3)
vector_de_co1_1=round(vector_de_co_1, digits = 3)

vector_de11=c(vector_de_co1,vector_de_pb1)
vector_de11_1=c(vector_de_co1_1,vector_de_pb1_1)
#View(vector_de11)


#Ahora crearemos los factores  vector_factores
#Creamos una columna para hacer la tabla de factores
v <- c("Control")
v_2 <- c("Control")
w <- rep(v, times = conteo_vector_de_co)
w_2 <- rep(v, times = conteo_vector_de_co_1)
#View(w)
#View(w_2)
ww=cbind(vector_de_co1,w)
#View(ww)
ww_2=cbind(vector_de_co1_1,w_2)
#View(ww_2)


v1 <- c("Pb")
v1_2 <- c("Pb")
w1 <- rep(v1, times = conteo_vector_de_pb)
w1_2 <- rep(v1, times = conteo_vector_de_pb_1)
#View(w1)



#View(ww1)

el_bind0_factores=c(w,w1)
#View(el_bind0_factores)
one=c(w_2,w1_2)
el_bind0_en=c(vector_de_co,vector_de_pb)
el_bind0_en_1=c(vector_de_co_1,vector_de_pb_1)
el_bind1_en=(el_bind0_en)
el_bind1_en_1=(el_bind0_en_1)
el_bind1_en_1_1=round(el_bind1_en, digits = 3)
two=round(el_bind1_en_1, digits = 3)
el_bind22=cbind(el_bind0_factores,el_bind1_en_1_1)
el_bind22_1=cbind(one,two)

#Calculo de Kruskall con datos sin log
# View(el_bind0_en)
# View(el_bind0_factores)
el_bind22_nvo=cbind(el_bind0_factores,el_bind0_en)
kruskal=kruskal.test(el_bind0_factores ~ el_bind0_en, data = el_bind22_nvo)
# Kruskal-Wallis chi-squared = 128549, df = 127210, p-value = 0.004064

el_bind2=as.data.frame(el_bind22)
el_bind2_1=as.data.frame(el_bind22_1)
#View(el_bind2)
#View(el_bind2_1)
write.csv(el_bind2_1, file="el_bind2_1.csv", quote = T, row.names = F)
exe10=("sed '1d' el_bind2_1.csv > el_bind2_1_2.csv")
system(exe10)
three=read.csv("el_bind2_1_2.csv", header = F, sep = ",")
#View(three)


#Grafico ponderado con expresion diferencial

library(ggplot2)
tiff("~/Escritorio/Mauricio/PB_SE/PB_SE/DE_PB.tiff",
     width=1200,height=820,units="px",
     pointsize=12,bg="white",res=300)  


ggplot(el_bind2, aes(x = el_bind0_factores, y = el_bind1_en)) +
  
  geom_boxplot()+
  #xlab('Treatments') +
  #dev.new(width=3, height=3)+
  xlab('') +
  ylab('log (N° ATP molecules)') +
  #Titulo del grafico
  #ggtitle('Code-energetic response to Pb') +
  theme_bw()+
  #tamaño del grafico.Proporcion
  coord_fixed(ratio = 0.3)+
  scale_y_log10()+

  
  theme (axis.text.x = element_text(colour="black", size=rel(1.2)))+
  theme(text = element_text(size=8))+
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.ticks.x = element_blank())



dev.off()

#Calculando diferencias
#KRUSKAL WALLIS
#Solo puedo calcularlas para replicas.
#kruskal.test(vector_de11 ~ vector_factores22, data = el_bind222)

sink("kruskal.txt") 
kruskal=kruskal.test(el_bind0_factores ~ el_bind1_en, data = el_bind2)
print(kruskal)
#Kruskal-Wallis rank sum test

#data:  el_bind0_factores by el_bind1_en
#Kruskal-Wallis chi-squared = 128549, df = 127848, p-value = 0.08301


tiff("~/Escritorio/Mauricio/PB_SE/PB_SE/unitario_PB.tiff",
     width=1200,height=820,units="px",
     pointsize=12,bg="white",res=300) 
####ALERTA !!!!!!  Ojo acà que puede venir duplicado el log
ggplot(three, aes(x = V1, y = V2)) +
  geom_bar(stat = "identity",width=0.5, position="dodge", fill="grey")+
  #xlab('Treatments') +
  #dev.new(width=3, height=3)+
  xlab('') +
  ylab('log (N° ATP molecules)') +
  #Titulo del grafico
  #ggtitle('Code-energetic response to Pb') +
  theme_bw()+
  #tamaño del grafico.Proporcion
  # coord_fixed(ratio = 4)+
  scale_y_log10()+
  
  
  theme (axis.text.x = element_text(colour="black", size=rel(1.2)))+
  theme(text = element_text(size=8))+
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.ticks.x = element_blank())
  
  # geom_boxplot()+
  # xlab('Trearments') +
  # ylab('log (N° ATP molecules)') +
  # ggtitle('Code-energetic response to Pb') +
  # theme_bw()+
  # theme(legend.position = "none",
  #       panel.grid = element_blank(),
  #       axis.ticks.x = element_blank())

dev.off()





# DIAGRAMA DE VENN
library(ggplot2)
library (VennDiagram)

exe2=("awk '{print $1}' control_sin_na_uno_log_energia_de.txt > nombres_solo_control")
system(exe2)
name_solo_cont=read.csv("nombres_solo_control", header = T,sep = "\t")
#View(name_solo_cont)
name_solo_cont_1=nrow(name_solo_cont)
# View(name_solo_cont_1)


exe3=("awk '{print $1}' pb_sin_na_uno_log_energia_de.txt > nombres_solo_pb")
system(exe3)
name_solo_pb=read.csv("nombres_solo_pb", header = T,sep = "\t")
name_solo_pb_1=nrow(name_solo_pb)
# View(name_solo_pb_1)


exe100 =("grep -wf nombres_solo_control nombres_solo_pb > grep_nombres_venn.csv")
system(exe100)
grep_2<-read.csv("grep_nombres_venn.csv", header = F,sep = "\t")
# View(grep_2)
grep_2_2=nrow(grep_2)
# View(grep_2_2)


png("~/Escritorio/Mauricio/PB_SE/PB_SE/VENN_Pb.png",
     width=1200,height=820,units="px",
     pointsize=12,bg="white",res=300)  

  
      draw.pairwise.venn (area1 = name_solo_pb_1,
                      area2 = name_solo_cont_1 ,
                      cross.area = grep_2_2,
                      #Nombres de los grupos
                      #  #Nombres de los grupos
                      category = c("Pb","Control"),
                      #  alpha = 0.75,
                      #  #Tamaño de numero
                      #  cex = 1.2,
                      #  #Tamaño letra
                      #  cat.cex = 1,
                      #  #Posicion de los nombres de los conjuntos
                       cat.pos = c(220, 140),
                      #  #Distancia de los nombres de los grupos al circulo
                      cat.dist = 0.1,
                      fontfamily=2,
                      #  #Grosor letra
                      cat.fontface = 1,
                      #  #Grosor numero
                      fontface = 1,
                      #  #Ancho linea conjuntos
                      #  sep.dist = 0.1,
                      lwd = 2,
                      #  #Separación de circulos
                      #  scaled = F
                   
      )
                                        
dev.off()

exe2=("awk '{print $1}' control_sin_na_uno_log_energia_de.txt > nombres_solo_control")
system(exe2)
name_solo_cont=read.csv("nombres_solo_control", header = T,sep = "\t")
#View(name_solo_cont)
name_solo_cont_1=nrow(name_solo_cont)
# View(name_solo_cont_1)


exe3=("awk '{print $1}' pb_sin_na_uno_log_energia_de.txt > nombres_solo_pb")
system(exe3)
name_solo_pb=read.csv("nombres_solo_pb", header = T,sep = "\t")
name_solo_pb_1=nrow(name_solo_pb)
# View(name_solo_pb_1)


exe100 =("grep -wf nombres_solo_control nombres_solo_pb > grep_nombres_venn.csv")
system(exe100)
grep_2<-read.csv("grep_nombres_venn.csv", header = F,sep = "\t")
# View(grep_2)
grep_2_2=nrow(grep_2)
# View(grep_2_2)












######
#ESTADISTICA PARA GRAFICA DE VENN
#####

#hice en la terminal

exe_a=("grep -wf nombres_solo_control nombres_solo_pb > compartidos_co_ca.csv")
system(exe_a)
grep_1=read.csv("compartidos_co_ca.csv", col.names = F)

compartidos=nrow(grep_1)

soloco1=read.csv("nombres_solo_control", col.names = F)
soloco=nrow(soloco1)


solopb1=read.csv("nombres_solo_pb", col.names = F)
solopb=nrow(solopb1)


solo_control= soloco - compartidos
solo_control
solo_plomo=solopb - compartidos
solo_plomo

#Solo control=990
#Solo plomo= 442

# Calculo para estimar la pbb de que los resultados anteriores se den por azar.

library("gmp")
#Sino se puede instalar
#sudo apt-get install libgmp-dev
#luego
#install.packages("gmp")

####Correr por default
#####Crear la función que calcula p = Σ (m,i)(N-m,n-i)/(N,n)
enrich_pvalue <- function(N, A, B, k)
{
  m <- A + k
  n <- B + k
  i <- k:min(m,n)
  
  as.numeric( sum(chooseZ(m,i)*chooseZ(N-m,n-i))/chooseZ(N,n) )
}

#Ingresar los valores
#N=14800, A=Total gupo a (menos los comunes), B = Total gupo B (menoslos comunes)
#K=comunes

#From 15220 genes, set A is 1850+195 genes,
#set B is 195+596 genes, overlap is 195 genes. Their p value is 2e-26.
#enrich_pvalue(15220, 1850, 596, 195)
sink("valorp_venn.txt") 
valorp=enrich_pvalue(29127, solo_control, solo_plomo, compartidos)
print(valorp)
#todos los genes involucrados  # genes solo co #genes solo trat #genes comunes

#Para los datos de Pb
#valorp=enrich_pvalue(29127, 990, 442, 20112)
#From 15220 genes, set A is 1850+195 genes, set B is 195+596 genes, overlap is 195 genes. Their p value is 2e-26.
#valorp1=enrich_pvalue(N=15220, A=1850, B=596, k=195)
#valorp1=enrich_pvalue(N=29127, A=990, B=442, k=20112)
#Resultado = 0

#summary(valorp)
valorp

#RESULTADO: 1.91221e-18

  
  
  #Numero esperado por azar
  sink("numero_esperado_por_azar.txt") 
  NEPA= solo_control*(solo_plomo/29127)
  print(NEPA)
  
  #Solo co(solo_cd/total)=
  #15
  
  
  
  
  sink("IC_CO_DE.txt") 
  b=wilcox.test(log(vector_de_co),
                alternative="two.sided",
                correct=TRUE,
                conf.int=TRUE,
                conf.level=0.95) 
  print(b)
  #Wilcoxon signed rank test with continuity correction
  
  
  
  #Intervalo de confianza para Pb-DE
  sink("IC_PB_DE.txt") 
  a<-wilcox.test(log(vector_de_pb),
                 alternative="two.sided",
                 correct=TRUE,
                 conf.int=TRUE,
                 conf.level=0.95) 
  
  print(a)
  
  #Intervalo de confianza para Co-Code
  sink("IC_CO_CODE.txt") 
  c=wilcox.test(log(vector_de_co_1),
                alternative="two.sided",
                correct=TRUE,
                conf.int=TRUE,
                conf.level=0.95)
  print(c)
  
  #Intervalo de confianza para Co-Code
  sink("IC_PB_CODE.txt") 
  d=wilcox.test(log(vector_de_pb_1),
                alternative="two.sided",
                correct=TRUE,
                conf.int=TRUE,
                conf.level=0.95) 
  print(d)