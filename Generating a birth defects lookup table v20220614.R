require(tidyverse)

# Variable names ----------------------------------------------------------

diagnoses <- c(
'birth.defect',

'cns.anomaly',
'anencephalus',
'spina.bifida.without.anencephalus',
'hydrocephalus.without.spina.bifida',
'encephalocele',
'microcephalus', 
'reduction.deformities.of.brain',

'eye.anomaly',
'anophthalmia.microphthalmia',
'congenital.cataract',
'aniridia', 

'ear.face.neck.anomaly',
'anotia.microtia.or.anom.causing.hearing.impairment',

'heart.circulatory.anom',
'common.truncus',
'tetralogy.of.fallot',
'tgv.dorv', # This includes transposition of the great vessels and double outlet right ventricle
'atrioventricular.septal.defect',
'total.anomalous.pulmonary.venous.return',
'hypoplastic.left.heart.syndrome',
'interrupted.aortic.arch',
'coarctation.of.aorta',
'aortic.valve.stenosis',
'pulm.valve.atresia.stenosis',
'ebstein.anomaly',
'pulm.artery.anomalies',
'ventricular.septal.defect',
'asd.or.pfo',
'single.ventricle',
'tricuspid.valve.atresia.stenosis',
'patent.ductus.arteriosus',

'respiratory.anomaly',
'lung.agenesis.hypoplasia',
'choanal.atresia',

'oral.clefts',
'cleft.palate.alone',
'cleft.lip.w.or.wo.cleft.palate',

'digestive.anomaly',
'esophageal.atresia.te.fistula',
'rectal.large.instestinal.atresia.stenosis',
'pyloric.stenosis',
'hirschsprung.disease',
'biliary.atresia',
'small.intestinal.atresia.stenosis',

'genitourinary,anomaly',
'renal.agenesis.or.hypoplasia',
'bladder.exstrophy',
'obstructive.genitourinary.defects',
'hypospadias',
'epispadias',
'cloacal.exstrophy',

'musculoskeletal.anomaly',
'limb.reduction.deformities',
'upper.limb.reduction.deformities',
'lower.limb.reduction.deformities',
'gastroschisis',
'omphalocele',
'congenital.hip.dislocation',
'diaphragmatic.hernia',
'clubfoot',
'craniosynostosis',

'integument.anomaly',

'genetic.anomaly',
'any.single.gene.anomaly',
'neurofibromatosis',
'tuberous.sclerosis',
'any.chromosomal.anomaly',
'trisomy.13',
'trisomy.18',
'trisomy.21',
'gonadal.dysgenesis',
'other.autosomal.deletions',
'deletion.13q',
'deletion.22q'
)

# BPA codes for each ------------------------------------------------------

bpa.codes <- c( 
str_c('^216','^237.7', '^7[45]', sep = '|'), # any defect
'^74[012]', # any CNS
'^740.[01]', # anencephalus
'^741', # spina bifida
'^742.3', # hydrocephalus
'^742.0[0-9]', # encephalocele
'^742.1', # microcephalus
'^742.2', # holoprosencephaly

'^743', # any eye 
'^743.[01]', # anophthalmia/microphthalmia
'^743.32', # congenital cataract
'^743.42', # aniridia

'^744.[0-9][0-9][012345679]', # any ear, face, or neck
str_c('^744.0', '^744.[02]1', sep = '|'), # anotia or microtia

'^74[567]', # any cardiovascular anomaly
'^745.0', # common truncus
str_c('^745.2', '747.31', sep = '|'), # Tetralogy of Fallot
'^745.1', # transposition of great vessels, including double outlet right ventricle
'^745.6', # AVSD (AKA endocardial cushion defect)
'^747.42', # TAPVR (total anomalous pulmonary venous return)
'^746.7', # hypoplastic left heart syndrome (HLHS)
'^747.21[567]', # interrupted aortic arch
'^747.1', # coarctation of the aorta
'^746.3', # aortic valve stenosis
'^746.0[01]', # pulmonary valve atresia/stenosis
'^746.2', # Ebstein anomaly
'^747.3', # pulmonary artery anomalies
str_c('^745.4[012]', '^745.48[01234569]', '^745.49[012345679]', sep = '|'), # ventricular septal defect
'^745.5', # atrial septal defect or patent foramen ovale
'^745.3', # single ventricle
'^746.10[^5]', # tricuspid valve atresia/stenosis
'^747.0', # patent ductus arteriosus

'^748', # any respiratory system anomaly
'^748.5', # lung agenesis
'^748.0', # choanal atresia

'^749', # orofacial clefts
'^749.0', # cleft palate without cleft lip
'^749.[12]', # cleft lip with or without cleft palate

'^75[01]', # digestive system anomalies
'^750.3[^8]', # esophageal atresia or tracheoesophageal fistula
'^751.2', # rectal or large intestinal atresia or stenosis
'^750.5', # pyloric stenosis
'^751.3', # Hirschsprung disease
'^751.65', # biliary atresia
'751.1', # small intestinal atresia or stenosis

'^75[23]', # any genitourinary anomaly
'^753.0', # renal agenesis or hypoplasia
'^753.5', # bladder exstrophy
'^753.[26]', # obstructive genitourinary defects
str_c('^752.60','^752.620', sep = '|'), # hypospadias
'^752.61', # epispadias
'^751.55', # cloacal exstrophy

'^75[456]', # conganomalies.musculoskeletal.sys
'^755.[234][0-9][^8]', # limb.reduction.deformities
'^755.2', # upper.limb.reduction.deformities
'^755.3', # lower.limb.reduction.deformities
'^756.71', # gastroschisis
'^756.70', # omphalocele
'^754.3', # congenital.hip.dislocation
'^756.61', # diaphragmatic.hernia
str_c('^754.50', '754.73', sep = '|'), # clubfoot
'^756.0[0-3]', # craniosynostosis

str_c('^757', '^216', '^228.0[01]', sep = '|'), # any integument anomaly

str_c('^237.7', '^758', '^759.[56]', sep = '|'), # any genetic anomaly 
str_c('^237.7','^759.[56]', sep = '|'), # any single gene anomaly (really NF or TSC)
'^237.70', # neurofibromatosis
'^759.5', # tuberous sclerosis complex
'^758', # any chromosomal anomaly
'^758.1', # trisomy 13 AKA Patau syndrome
'^758.2', # trisomy 18 AKA Edwards syndrome
'^758.0', # trisomy 21
'^758.6', # gonadal dysgenesis (mostly Turner syndrome)
'758.3', # other autosomal deletion syndromes
'^758.33', # deletion of 13q
'^758.37' # deletion of 22q
)

# "Canonical" BPA codes from registry docs ------------------------------------

#' This section provides more precise ("canonical") code ranges that would be expected from registries using the CDC-BPA system.
#' It's deprecated, largely due to constraints on my time and the trouble of maintaining it.

#canonical.bpa.codes <- c(
#'^74[012].[0-9][0-9][^8]', # any CNS
#paste('^740.0[01238][0-9]', '^740.10[0-9]', sep = ', '), # anencephalus
#paste('^741.0[0-7][0-9]', '^741.08[0567]', '^741.09[0-9]', sep = ', '), # spina bifida
#'^742.3[01289][0-9][^8]', # hydrocephalus
#paste('^742.0[09][0-9]', '^742.08[056][^8]', sep = ', '), # encephalocele
#'^742.10[0-9]', # microcephalus
#'^742.26[0-9]0', # holoprosencephaly

#'^743.[0-9][0-9][^8]', # any eye
#'^743.[01]0[^8]', # anophthalmia/microphthalmia
#'^743.32[056]', # congenital cataract
#'^743.42[^8]', # aniridia

#'^744.[0-9][0-9][^8]', # any ear, face, or neck
#paste('^744.0[01239][^8]', '^744.[02]1[012345679]', sep = ','), # anotia, microtia, or anomalies causing hearing impairment

#'^74[567].[0-9][0-9][^8]', # any cardiovascular
#'^745.0[01][012345679]', # common truncus
#paste('^745.20[012345679]', '747.31[012345679]' , sep =','), # Tetralogy of Fallot                     
#'^745.1[01289][^8]', #transposition of great vessels
#'^745.6[01389][^8]', # AVSD (AKA endocardial cushion defect)
#'^747.42[^8]', # TAPVR (total anomalous pulmonary venous return)
#'^746.7[0][^8]', # hypoplastic left heart syndrome (HLHS)
#'^747.21[567]', # interrupted aortic arch
#'^747.1[019][^8]', # coarctation of the aorta
#'^746.30[^8]', # aortic valve stenosis
#'^746.0[01][^8]', # pulmonary valve atresia/stenosis,
#'^746.20[^8]', # Ebstein anomaly
#paste('^747.3[0123489]0', '^747.325', sep = ','), # pulmonary artery anomalies
#paste('^745.4[012][^8]', '^745.48[01234569]', '^745.49[012345679]', sep = ','), # ventricular septal defect
#'^745.5[01289]', # atrial septal defect or patent foramen ovale
#'^745.30[^8]', # single ventricle
#'^746.10[06]', # tricuspid valve atresia/stenosis
#'^747.00[^8]', # patent ductus arteriosus

#'^748.[0-9][0-9][^8]', # any respiratory system anomaly
#'^748.5[01289][^8]', # lung agenesis
#'^748.00[^8]', # choanal atresia

#'^749.[012][0-9][^8]', # orofacial clefts
#'^749.0[02346789][^8]', # cleft palate without cleft lip
#'^749.[12][02]', # cleft lip with or without cleft palate

#'^75[01].[0-9][0-9][^8]', # digestive system anomalies
#paste('^750.3[01345]0', '^750.32[05]', sep = ','), # esophageal atresia or tracheoesophageal fistula
#'^751.2[0-4][^8]', # rectal or large intestinal atresia or stenosis
#'^750.5[018][^8]', # pyloric stenosis
#'^751.3[0-4][^8]', # Hirschsprung disease
#'^751.65[^8]', # biliary atresia
#paste('751.1[012][^8]', '^751.19[05][^8]', sep = ','), # small intestinal atresia or stenosis

#'^75[23].[0-9][0-9][^8]', # any genitourinary anomaly
#'^753.0[01][09]', # renal agenesis or hypoplasia
#'^753.50[^8]', # bladder exstrophy
#paste('^753.2[0129][^8]', '^753.6[01239][^8]', sep = ',' ), # obstructive genitourinary defects
#paste('^752.60[0567]','^752.620', sep = ','), # hypospadias
#'^752.61[^8]', # epispadias
#'^751.550', # cloacal exstrophy

#'^75[456].[0-9][0-9][^8]', # conganomalies.musculoskeletal.sys',
#paste('^755.2[01234579][^8]', '^755.2[68][05]', '^755.3[0123459][^8]', '^755.3[68][056]', '^755.4[0123489][^8]', sep = ',' ), # limb.reduction.deformities',
#paste('^755.2[01234579][^8]', '^755.2[68][05]', sep = ','), # upper.limb.reduction.deformities',
#paste('^755.3[0123459][^8]', '^755.3[68][056]', sep = ','), # lower.limb.reduction.deformities',
#'^756.71[^8]', # gastroschisis',
#'^756.70[^8]', # omphalocele',
#'^754.3[01][^8]', # congenital.hip.dislocation
#'^756.61[^8]', # diaphragmatic.hernia
#paste('^754.50[^8]', '754.73[^8]', sep = ','), # clubfoot
#'^756.0[0123]', # craniosynostosis

#paste('^757.[0-9][0-9][^8]', '^216.[0-9][0-9][^8]', '^228.0[01][^8]', sep = ','), # any integument anomaly

#paste('^237.7', '^758', '^759.50', '^759.6[012389]', sep = ','), # any genetic anomaly 
#paste('^237.70', '^759.[56]', sep = ','), # any single gene anomaly
#'^237.70', # neurofibromatosis,
#'^759.50', # tuberous sclerosis complex
#'^758.[0-9][0-9][^8]', # any chromosomal anomaly
#'^758.1[012349]', # trisomy 13 AKA Patau syndrome
#paste('^758.2[01234]', '^758.29[05]', sep = ','), # trisomy 18 AKA Edwards syndrome
#'^758.0[012349]', # trisomy 21
#'^758.6[019]', # gonadal dysgenesis (mostly Turner syndrome)
#'^758.33', # deletion of 13q
#'^758.37' # deletion of 22
#)

# ICD9 codes --------------------------------------------------------------

icd9.codes <- c( 
  str_c('^216', '^237.7', '^7[45]', sep = '|'), # any birth defect
  '^74[012]', # any CNS
  '^740.[01]', # anencephalus
  '^741', # spina bifida
  '^742.3', # hydrocephalus
  '^742.0', # encephalocele
  '^742.1', # microcephalus
  '^742.2', # reduction deformities of the brain, primarily holoprosencephaly
  
  '^743', # any eye 
  '^743.[01]', # anophthalmia/microphthalmia
  '^743.3[0-4]', # congenital cataract
  '^743.45', # aniridia
  
  '^744', # any ear, face, or neck
  str_c('^744.01', '^744.23', sep = "|"), # anotia or microtia
  
  '^74[567]', # any cardiovascular anomaly
  '^745.0', # common truncus
  '^745.2', # Tetralogy of Fallot
  '^745.1', # transposition of great vessels, including double outlet right ventricle
  '^745.6', # AVSD (AKA endocardial cushion defect)
  '^747.41', # TAPVR (total anomalous pulmonary venous return)
  '^746.7', # hypoplastic left heart syndrome (HLHS)
  '^747.11', # interrupted aortic arch
  '^747.10', # coarctation of the aorta
  '^746.3', # aortic valve stenosis
  '^746.0[12]', # pulmonary valve atresia/stenosis
  '^746.2', # Ebstein anomaly
  '^747.3', # pulmonary artery anomalies
  '^745.4', # ventricular septal defect
  '^745.5', # atrial septal defect or patent foramen ovale
  '^745.3', # single ventricle
  '^746.1', # tricuspid valve atresia/stenosis
  '^747.0', # patent ductus arteriosus
  
  '^748', # any respiratory system anomaly
  '^748.5', # lung agenesis
  '^748.0', # choanal atresia
  
  '^749', # orofacial clefts
  '^749.0', # cleft palate without cleft lip
  '^749.[12]', # cleft lip with or without cleft palate
  
  '^75[01]', # digestive system anomalies
  '^750.3', # esophageal atresia or tracheoesophageal fistula
  '^751.2', # rectal or large intestinal atresia or stenosis
  '^750.5', # pyloric stenosis
  '^751.3', # Hirschsprung disease
  '^751.61', # biliary atresia
  '751.1', # small intestinal atresia or stenosis
  
  '^75[23]', # any genitourinary anomaly
  '^753.0', # renal agenesis or hypoplasia
  '^753.5', # bladder exstrophy
  '^753.[26]', # obstructive genitourinary defects
  '^752.61', # hypospadias
  '^752.62', # epispadias
  '^751.5', # cloacal exstrophy
  
  '^75[456]', # conganomalies.musculoskeletal.sys
  '^755.[234]', # limb.reduction.deformities
  '^755.2', # upper.limb.reduction.deformities
  '^755.3', # lower.limb.reduction.deformities
  '^756.73', # gastroschisis
  '^756.72', # omphalocele
  '^754.3[01]', # congenital.hip.dislocation
  '^756.6', # diaphragmatic.hernia
  str_c('^754.51', '754.70', sep = '|'), # clubfoot
  '^756.0', # craniosynostosis
  
  str_c('^757', '^216', '^228.0', sep = '|'), # any integument anomaly
  
  str_c('^237.7', '^758', '^759.5', sep = '|'), # any genetic anomaly 
  str_c('^237.7','^759.5', sep = '|'), # any single gene anomaly (really NF or TSC)
  '^237.7', # neurofibromatosis
  '^759.5', # tuberous sclerosis complex
  '^758', # any chromosomal anomaly
  '^758.1', # trisomy 13 AKA Patau syndrome
  '^758.2', # trisomy 18 AKA Edwards syndrome
  '^758.0', # trisomy 21
  '^758.6', # gonadal dysgenesis (mostly Turner syndrome)
  '^758.3', # other autosomal deletion syndromes
  '^999.999', # deletion of 13q
  '^999.999' # deletion of 22q
)

# ICD10 codes -------------------------------------------------------------

icd10.codes <- c( 
  '^Q[0-9]', # birth defect
  
  '^Q0[0-7]', # any CNS
  '^Q00.[01]', # anencephalus
  '^Q05', # spina bifida
  '^Q03', # hydrocephalus
  '^Q01', # encephalocele
  '^Q02', # microcephalus
  '^Q04.[0123]', # reduction deformities of the brain, primarily holoprosencephaly
  
  '^Q1[0-5]', # any eye 
  '^Q11.[0-3]', # anophthalmia/microphthalmia
  '^Q12.0', # congenital cataract
  '^Q13.1', # aniridia
  
  '^Q1[6-8]', # any ear, face, or neck
  '^Q17.[12]', # anotia or microtia
  
  '^Q2[0-8]', # any cardiovascular anomaly
  '^Q20.0', # common truncus
  '^Q21.3', # Tetralogy of Fallot
  '^Q20.3', # transposition of great vessels, including double outlet right ventricle
  '^Q21.2', # AVSD (AKA endocardial cushion defect)
  '^Q26.2', # TAPVR (total anomalous pulmonary venous return)
  '^Q23.4', # hypoplastic left heart syndrome (HLHS)
  '^Q25.21', # interrupted aortic arch
  '^Q25.1', # coarctation of the aorta
  '^Q23.0', # aortic valve stenosis
  '^Q22.[01]', # pulmonary valve atresia/stenosis
  '^Q22.5', # Ebstein anomaly
  str_c('^Q25.[56]', '^Q25.7[1-2]', '^Q25.79', sep = '|'), # pulmonary artery anomalies
  '^Q21.0', # ventricular septal defect
  '^Q21.1', # atrial septal defect or patent foramen ovale
  '^Q20.4', # single ventricle
  '^Q22.4', # tricuspid valve atresia/stenosis
  '^Q25.0', # patent ductus arteriosus
  
  '^Q3[0-4]', # any respiratory system anomaly
  '^Q33.[36]', # lung agenesis
  '^Q30.0', # choanal atresia
  
  '^Q3[5-7]', # orofacial clefts
  '^Q35', # cleft palate without cleft lip
  str_c('^Q36.[0-1]', '^Q36.9', '^Q37.[0-5]','^Q37.[8-9]', sep = '|'), # cleft lip with or without cleft palate
  
  str_c('^Q3[8-9]', '^Q4[0-5]', sep = '|'), # digestive system anomalies
  '^Q39.[0-2]', # esophageal atresia or tracheoesophageal fistula
  '^Q42', # rectal or large intestinal atresia or stenosis
  '^Q40.0', # pyloric stenosis
  '^Q43.1', # Hirschsprung disease
  '^Q44.2', # biliary atresia
  '^Q41', # small intestinal atresia or stenosis
  
  str_c('^5[0-6]', '^Q6[0-4]', sep = '|'), # any genitourinary anomaly
  '^Q60', # renal agenesis or hypoplasia
  '^Q64.1[29]', # bladder exstrophy
  str_c('^Q62.0', '^Q62.1[0-2]', '^Q62.2', '^Q62.3[1-2]', '^Q64.2', '^Q64.3[1-3]', '^Q64.39', sep = '|'), # obstructive genitourinary defects
  str_c('^Q54.[0-3]', '^Q54.[8-9]', sep = '|'), # hypospadias
  '^Q64.0', # epispadias
  '^Q43.[5-7]', # cloacal exstrophy
  
  str_c('^6[5-9]', '^Q7[0-4]', sep = '|'), # conganomalies.musculoskeletal.sys
  '^Q7[1-3]', # limb.reduction.deformities
  '^Q71', # upper.limb.reduction.deformities
  '^Q72', # lower.limb.reduction.deformities
  '^Q79.3', # gastroschisis
  '^Q79.2', # omphalocele
  '^Q65.[0-5]', # congenital.hip.dislocation
  '^Q79.0', # diaphragmatic.hernia
  '^Q66.0[0-2]', # clubfoot
  '^Q75.0', # craniosynostosis
  
  '^Q8[0-4]', # any integument anomaly
  
  str_c('^Q85.0','^Q85.1', '^Q91.[0-7]', '^Q90.[0129]', '^Q96', '^Q93', sep = '|'), # any genetic anomaly 
  str_c('^Q85.0','^Q85.1', sep = '|'), # any single gene anomaly (really NF or TSC)
  '^Q85.0', # neurofibromatosis
  '^Q85.1', # tuberous sclerosis complex
  str_c('^Q91.[0-7]', '^Q90.[0129]', '^Q96', '^Q93', sep = ','), # any chromosomal anomaly
  '^Q91.[4-7]', # trisomy 13 AKA Patau syndrome
  '^Q91.[0-3]', # trisomy 18 AKA Edwards syndrome
  '^Q90.[0129]', # trisomy 21
  '^Q96', # gonadal dysgenesis (mostly Turner syndrome)
  '^Q93', # other autosomal deletion syndromes
  '^Q93.89', # deletion of 13q
  '^Q93.81' # deletion of 22q
)

# The organ systems to which individual birth defects belong --------------

#' NA for organ system level variables. Will be used later within the compute.bd.vars function.
systems <- c(
rep(NA, 2),
rep('cns', 6),
NA,
rep('eye', 3),
NA,
'ear.face.neck',
NA,
rep('heart.circ.sys', 17),
NA,
rep('resp.sys', 2),
NA,
rep('clefts', 2),
NA,
rep('digestive.sys', 6),
NA,
rep('genitourinary', 6),
NA,
rep('musculoskeletal.sys', 9),
rep(NA, 2),
rep('syndromes',11)
)

# Combine each into a data frame ------------------------------------------

defects <- data.frame(diagnosis = diagnoses, 
                      bpa.codes = bpa.codes,
                      icd9.codes = icd9.codes,
                      icd10.codes = icd10.codes,
                      organ.system = systems,
                      stringsAsFactors = F)

saveRDS(defects, 
     file = '//smb-main.ad.bcm.edu/genepi2/Old_genepi2/Jeremy/GOBACK/Datasets/Expanded datasets/birth.defects.lookup.table.v20220614.rds')

write_csv(defects,
          file = '//smb-main.ad.bcm.edu/genepi2/Old_genepi2/Jeremy/GOBACK/Datasets/Expanded datasets/birth.defects.lookup.table.v20220614.csv')

rm(list = ls()); gc()
