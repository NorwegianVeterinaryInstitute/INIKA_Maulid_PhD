
library(tidyverse)

test_CPE<-joined_data%>%
  rename(CPE ='COLONY MORPHOLOGY ON CARBA.x')%>%
  filter(CPE == "Pinkish/Reddish")
# also three that originally are K.pneumoniae ( I guess thes ehave been selected)
# However the joined dataset shows some dscrepancies and therefore we need to gonthrough the cleaning and joineing steps once more!

# Note there are 11 CPE E.coli isolates, but one of these has been identified also on the CRGR3 plate,
# which makes it difficult to identify if this is another isoalte or not has only one Isolates_ID!

testECO<-joined_data%>%
filter(PROTOCOL == "CGR3")%>%
  rename(CRGR ='COLONY MORPHOLOGY ON C3GR.x')%>%
  select(TVLA_ID, Isolate.x,VITEK_MS_Results,CRGR, PROTOCOL, CTX_ED5, AMX_ED10)%>%
  filter(VITEK_MS_Results=="Escherichia coli")%>%
  
# BotH there is one isolate that has a mm zone of 30 of these isolates above so to be strict only 128 isolates should be reported as ESCR resistant E.coli

# 21 of the isolates that originally were K.pneumoniae have been conformed by VITEK, 
#  we decided to not include these inthe AST table as these might have been contaminated. ( or not pure samples)- 
#  Some of these ( of the selected ones) turned out to be K.pneumoinae by Malditof at NVI!


testKLEB<-joined_data%>%
  filter(PROTOCOL == "CGR3")%>%
  rename(CRGR ='COLONY MORPHOLOGY ON C3GR.x')%>%
  select(TVLA_ID, Isolate.x,VITEK_MS_Results, PROTOCOL, CTX_ED5, AMX_ED10)%>%
  filter(VITEK_MS_Results=="Klebsiella pneumoniae")

# 50 isolates identified by both the biochemical tests and confirmed by VITEK MS
#  8 isolates confirmed by VITEK as K.pneumoniae but biochemical tests identified these as E.coli, thus we decided to not include these isoaltes in the AST results

testSalmonella<- joined_data%>%
  filter(PROTOCOL == "SALM")%>%
  select(TVLA_ID, Isolate.x,VITEK_MS_Results, PROTOCOL, CTX_ED5, AMX_ED10, CITRATE.x, `TSI Media_Slope.x`,
         `TSI Media_Butt.x`, `TSI Media_Gas.x`, TSIMedia_H2S.x, UREASE.x)

# You need to make a case_when for to correct the results in th eoriginal file before the joining
# Citrate and Urease when it is  "d" it should be + for Citrate and for Urease it should be -, need 2 separate case_when for this.

  

filter(CRG3 == "Pinkish/Reddish")



test<-joined_data%>%
  filter( )%>%
  select(TVLA_ID, Isolate.x,VITEK_MS_Results, PROTOCOL)
mutate(CRG3 = 'COLONY MORPHOLOGY ON C3GR.x')%>%
  filter(CRG3 == "Pinkish/Reddish")
  
str(joined_data)