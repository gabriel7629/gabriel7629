chemin_acces<-"//calebasse/mariett231/Documents/stat chauvel/TP cœurs Poulet-20240419"

data<-read.table(paste(chemin_acces,"Data_fictives.txt",sep=""),sep="\t",header=T)

tdata<-read.table(paste(chemin_acces,"Data_fictives_t.txt",sep=""),sep="\t",header=T)

###Statistiques descriptives ---------------------------------
tdata$Traitement<-factor(tdata$Traitement,levels=c("T_AChE","AChE_1e-6","AChE_1e-5","AChE_1e-4","T_NORa","NORa_1e-6","NORa_1e-5","NORa_1e-4","T_X","X_concC" )) # Pour mettre les données dans l'ordre souhaité

par(mar=c(6, 4, 4, 2) + 0.1) # Pour augmenter les marges de la fenêtre graphique
boxplot(tdata$FC~tdata$Traitement,las=3,xlab="",ylab="FC (bpm)",col=c("#FFFFFF", "#AFA1FF", "#674CFF", "#2200FF", "#FFFFFF", "#FA9E9E", "#FF4F4F", "#FF0000", "#FFFFFF", "#E900FA"),font.title=2)
dev.off()# Pour fermer le plot et remettre les paramètres graphiques à 0

### intervalles confiance ---------------------------------

mTACHE <- mean(data$T_AChE)
sTACHE <- sd(data$T_AChE)
mTNORA <- mean(data$T_NORa)
sTNORA <- sd(data$T_NORa)
mTX <- mean(data$T_X)
sTX <- sd(data$T_X)

m10.4ACHE <- mean(data$AChE_1e.4)
s10.4ACHE <- sd(data$AChE_1e.4)
m10.4NORA <- mean(data$NORa_1e.4)
s10.4NORA <- sd(data$NORa_1e.4)
mC.X <- mean(data$X_concC)
sC.X <- sd(data$X_concC)

t1a2<-qt(0.975,49)

# TACHE
mTACHE-t1a2*(sTACHE/sqrt(50)) #borne inférieure
mTACHE+t1a2*(sTACHE/sqrt(50)) #borne supérieure

# TNORA
mTNORA-t1a2*(sTNORA/sqrt(50)) #borne inférieure
mTNORA+t1a2*(sTNORA/sqrt(50)) #borne supérieure

# TX
mTX-t1a2*(sTX/sqrt(50)) #borne inférieure
mTX+t1a2*(sTX/sqrt(50)) #borne supérieure

# 10.4ACHE
m10.4ACHE-t1a2*(s10.4ACHE/sqrt(50)) #borne inférieure
m10.4ACHE+t1a2*(s10.4ACHE/sqrt(50)) #borne supérieure

# 10.4NORA
m10.4NORA-t1a2*(s10.4NORA/sqrt(50)) #borne inférieure
m10.4NORA+t1a2*(s10.4NORA/sqrt(50)) #borne supérieure

# C.X
mC.X-t1a2*(sC.X/sqrt(50)) #borne inférieure
mC.X+t1a2*(sC.X/sqrt(50)) #borne supérieure

### Statistiques inférentielles ---------------------------------

#O1
#Test t de Student
#test d'homogénéité pour échantillons indépendants
#test bilatéral
#H0 : il n'existe pas de différence significative des FC entre les deux témoins ; HA : il existe une différence significative des FC entre les deux témoins
#formule : abs(mTACHE-mTNORA)/sqrt(s²pd/nACHE+s²pd/nNORA)

s2pd<-(49*sTACHE^2+49*sTNORA^2)/(50+50-2)

tcalc<-abs(mTACHE-mTNORA)/sqrt(s2pd/50+s2pd/50)
tseuil<-qt(0.025,df=(50+50-2))
tcalc>tseuil
tcalc<(-tseuil) #ou
tseuil2<-qt(0.975,df=(50+50-2))

tcalc<tseuil2 #revient au même
tcalc>-tseuil2 #revient au même

#Les deux conditions nécéssaires à H0 sont vérifiées, donc H0 ne peut pas être rejetée, il n'y a pas de différence au seuil de rejet de 5% (alpha = 0.05), les FC des témoins ne sont donc pas significativement différentes

#Calcul de la p-value
2*pt(q=-abs(tcalc), df=(50+50-2))#en bilatéral, ce que l'on cherche ici. p >0,05, on ne rejette donc pas H0
#ou
2*pt(q=(tcalc), df=(50+50-2),lower.tail = F)

#Formule R
t.test(data$T_AChE,data$T_NORa,var.equal = T)

#O2
#Test t de Student
#test de conformité
#test bilatéral (la FC du témoin + 25% n'est pas censée être différente de celle mesurée pour X_ConcC)
#H0 : il n'existe pas de différence significative des FC entre témoin et coeurs en contact avec du X à C ; HA : il existe une différence significative des FC entre le témoin et coeurs en contact avec du X à C
#formule : (mTX-25%mTX)/(sTX/sqrt(50))

mTX25<-mTX*1.25 # calcul d'une FC témoin 25% plus rapide

tcalc<-(mC.X-mTX25)/(sC.X/sqrt(50))
tseuil<-qt(0.025,df=49) # Test bilatéral, les coeurs traités avec X_conC sont censés présenter une moyenne égale à mTx25
tcalc>tseuil
tcalc<(-tseuil)
#Donc H0 est acceptée, il n'y a pas de différence entre le témoin augmenté de 25% et les coeurs traités avec X concentrée à C

#Calcul de la p-value

2*pt(q=-abs(tcalc), df=(50-1))#p < 0.05, H0 est donc bien rejetée

# Formule R
t.test(data$X_concC,mu=mTX25,var.equal = T)

#O3

#T_AChE vs AChE.1.10-4 en bilatéral

#Test t de Student
#test d'homogénéité pour échantillons appariés
#test bilatéral
#H0 : il n'existe pas de différence significative des FC entre témoin et coeurs en contact avec AChE à 1.10-4 ; HA : il existe une différence significative des FC entre le témoin et coeurs en contact avec AChE à 1.10-4
#formule : mdiff/(sdiff/sqrt(50))

mdiff<-mean(data$T_AChE-data$AChE_1e.4)
sdiff<-sd(data$T_AChE-data$AChE_1e.4)

tcalc<-mdiff/(sdiff/sqrt(50))
tseuil<-qt(0.025,df=(49))
tcalc>tseuil
tcalc<(-tseuil)
#L'une des deux conditions n'est pas vérifiée, donc H0 est rejetée, il y a bien une différence significative entre la moyenne du témoin et celle des coeurs en contact avec l'AChE à 1.10-4

#Calcul de la p-value

2*pt(q=-abs(tcalc), df=(50-1),lower.tail=T)#p <<< 0.05, H0 est donc bien rejetée
#ou
2*pt(q=(tcalc), df=(50-1),lower.tail=F)

# Formule R
t.test(data$T_AChE,data$AChE_1e.4,paired = T,var.equal = T)

#T_AChE vs AChE.1.10-4 en unilatéral

#Test t de Student
#test d'homogénéité pour échantillons appariés
#test unilatéral à gauche (on suppose que l'acétylcholine réduit le nb de bpm)
#H0 : la moyenne du nombre de bpm des coeurs en contact avec AChE à 1e-4 n'est pas inférieure (peut-être égale ou supérieure) à celle des coeurs en condition témoin ; HA : la moyenne du nombre de bpm des coeurs en contact avec AChE à 1e-4 est inférieure à celle des coeurs témoins
#formule : mdiff/(sdiff/sqrt(50))

mdiff<-mean(data$AChE_1e.4-data$T_AChE)
sdiff<-sd(data$AChE_1e.4-data$T_AChE)

tcalc<-mdiff/(sdiff/sqrt(50))
tseuil<-qt(0.05,df=49)
tcalc>tseuil
#Donc H0 est rejetée, le nombre de bpm des coeurs en contact avec l'AChE à 1e-4 est inférieur à celui des coeurs en condition témoin (ralentissement par AChE)

#Calcul de la p-value

pt(q=tcalc, df=(50-1),lower.tail = T)#p <<< 0.05, H0 est donc bien rejetée

# Formule R
t.test(data$AChE_1e.4,data$T_AChE,paired = T,alternative="less",var.equal = T)

#T_NORa vs NORa.1.10-4 en bilatéral

#Test t de Student
#test d'homogénéité pour échantillons appariés
#test bilatéral
#H0 : il n'existe pas de différence significative des FC entre témoin et coeurs en contact avec NORa à 1.10-4 ; HA : il existe une différence significative des FC entre le témoin et coeurs en contact avec NORa à 1.10-4
#formule : mdiff/(sdiff/sqrt(50))

mdiff<-mean(data$T_NORa-data$NORa_1e.4)
sdiff<-sd(data$T_NORa-data$NORa_1e.4)

tcalc<-mdiff/(sdiff/sqrt(50))
tseuil<-qt(0.025,df=(50-1))
tcalc>tseuil
tcalc<(-tseuil)
#L'une des deux conditions n'est pas repsectée, donc H0 est rejetée, il y a bien une différence significative entre la moyenne du témoin et celle des coeurs en contact avec l'NORa à 1.10-4

#Calcul de la p-value

2*pt(q=-abs(tcalc), df=(50-1))#p <<< 0.05, H0 est donc bien rejetée

# Formule R
t.test(data$T_NORa,data$NORa_1e.4,paired = T,var.equal = T)


#T_NORa vs NORa.1.10-4 en uniléral

#Test t de Student
#test d'homogénéité pour échantillons appariés
#test unilatéral à droite (on suppose que la noradrénaline augmente le nb de bpm)
#H0 : le nombre de bpm des coeurs en contact avec NORa à 1e-4 ne sont pas supérieurs à ceux du témoins ; HA : le nombre de bpm des coeurs en contact avec NORa à 1e-4 sont supérieurs à ceux du témoins
#formule : abs(mdiff)/(sdiff/sqrt(50))

mdiff<-mean(data$NORa_1e.4-data$T_NORa)
sdiff<-sd(data$T_NORa-data$NORa_1e.4)

tcalc<-abs(mdiff)/(sdiff/sqrt(50))
tseuil<-qt(0.95,df=(49))
tcalc<tseuil

#Tests synonymes :
tseuil2<-qt(0.05,df=(49),lower.tail=F)
tcalc<tseuil2

tseuil3<-qt(0.05,df=(49),lower.tail=T)
tcalc<(-tseuil3)


#Donc H0 est rejetée, le nombre de bpm des coeurs en contact avec la NORa à 1e-4 est supérieur à celui des coeurs en condition témoin (stimulation pas NORa)

#Calcul de la p-value

pt(q=tcalc, df=(50-1),lower.tail=F)#p <<< 0.05, H0 est donc bien rejetée

# Formule R
t.test(data$NORa_1e.4,data$T_NORa,alternative="greater",paired = T,var.equal = T)


#Autres tests (détails non exigés)

#AchE
#en contexte bilatéral :

t.test(data$AChE_1e.5,data$T_AChE,paired = T,var.equal = T)
t.test(data$AChE_1e.6,data$T_AChE,paired = T,var.equal = T)

# // // unilatéral
t.test(data$AChE_1e.5,data$T_AChE,alternative="less",paired = T,var.equal = T)
t.test(data$AChE_1e.6,data$T_AChE,alternative="less",paired = T,var.equal = T)

#NorA
#en contexte bilatéral :

t.test(data$NORa_1e.5,data$T_NORa,paired = T,var.equal = T)
t.test(data$NORa_1e.6,data$T_NORa,paired = T,var.equal = T)

# // // unilatéral
t.test(data$NORa_1e.5,data$T_NORa,alternative="greater",paired = T,var.equal = T)
t.test(data$NORa_1e.6,data$T_NORa,alternative="greater",paired = T,var.equal = T)
