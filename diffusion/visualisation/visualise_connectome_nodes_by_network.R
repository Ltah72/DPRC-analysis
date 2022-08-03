


pacman::p_load(ggplot2, corrplot, viridis, magrittr)


setwd('V:/Vault/NECTAR_data/LENORE/derivatives/groups/F0/diff_data/cross-sectional/connectome/stats_results/weighted/ConnectomeWhole_linear_trend_stats')

#read in the family-wise error (fwe) stat files
connectome_whole_360_nodes <- read.csv('outputWhole_360Nodes_connectome_linear_trend_Rvisualise_connectome_fwe_1mpvalue_t1.csv', header = FALSE)
#remove the rows and columns if they contain all zeros in them.
#z<-y[rowSums(y[])>0,colSums(y[])>0]



#nodes <- c(73,253,67,247,97,277,98,278,26,206,70,250,71,251,87,267,68,248,83,263,85,265,84,264,86,266,40,220,41,221,55,235,44,224,43,223,36,216,39,219,37,217,48,228,95,275,49,229,117,297,50,230,47,227,42,222,45,225,46,226,29,209,143,323,151,331,150,330,149,329,148,328,116,296,147,327,146,326,145,325,144,324)
#nodes_sorted <- sort(nodes)

nodes_network_ordered<-c(1,	181,	4,	5,	6,	184,	185,	186,	3,	13,	16,	17,	19,	152,	183,	193,	196,	197,	199,	332,	7,	18,	22,	153,	154,	160,	163,	187,	198,	202,	333,	334,	340,	343,	2,	20,	21,	23,	138,	156,	157,	158,	159,	182,	200,	201,	203,	318,	336,	337,	338,	339,	8,	9,	51,	52,	53,	188,	189,	231,	232,	233,	36,	37,	38,	39,	40,	41,	43,	44,	55,	216,	217,	218,	219,	220,	221,	223,	224,	235,	10,	11,	12,	54,	56,	78,	96,	190,	191,	192,	234,	236,	258,	276,	99,	100,	101,	102,	113,	279,	280,	281,	282,	293,	24,	103,	104,	105,	124,	173,	174,	204,	283,	284,	285,	304,	353,	354,	107,	123,	125,	128,	129,	130,	175,	176,	287,	303,	305,	308,	309,	310,	355,	356,	106,	108,	109,	110,	111,	112,	114,	115,	167,	168,	169,	178,	286,	288,	289,	290,	291,	292,	294,	295,	347,	348,	349,	358,	118,	119,	120,	122,	126,	127,	135,	155,	298,	299,	300,	302,	306,	307,	315,	335,	131,	132,	133,	134,	136,	137,	172,	177,	311,	312,	313,	314,	316,	317,	352,	357,	25,	28,	139,	140,	141,	205,	208,	319,	320,	321,	29,	42,	45,	46,	47,	48,	49,	50,	95,	117,	209,	222,	225,	226,	227,	228,	229,	230,	275,	297,	116,	143,	144,	145,	146,	147,	148,	149,	150,	151,	296,	323,	324,	325,	326,	327,	328,	329,	330,	331,	14,	15,	27,	30,	31,	32,	33,	34,	35,	121,	142,	161,	162,	194,	195,	207,	210,	211,	212,	213,	214,	215,	301,	322,	341,	342,	57,	58,	59,	60,	61,	62,	63,	64,	65,	69,	88,	164,	165,	166,	179,	180,	237,	238,	239,	240,	241,	242,	243,	244,	245,	249,	268,	344,	345,	346,	359,	360,	66,	72,	89,	90,	91,	92,	93,	94,	170,	246,	252,	269,	270,	271,	272,	273,	274,	350,	74,	75,	76,	77,	79,	80,	81,	82,	171,	254,	255,	256,	257,	259,	260,	261,	262,	351,	26,	67,	68,	70,	71,	73,	83,	84,	85,	86,	87,	97,	98,	206,	247,	248,	250,	251,	253,	263,	264,	265,	266,	267,	277,	278)
nodes_sorted_numerical <- sort(nodes_network_ordered)

#section_22_cortices <- c(1,	1,	2,	2,	2,	2,	2,	2,	3,	3,	3,	3,	3,	3,	3,	3,	3,	3,	3,	3,	4,	4,	4,	4,	4,	4,	4,	4,	4,	4,	4,	4,	4,	4,	5,	5,	5,	5,	5,	5,	5,	5,	5,	5,	5,	5,	5,	5,	5,	5,	5,	5,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	7,	7,	7,	7,	7,	7,	7,	7,	7,	7,	7,	7,	7,	7,	7,	7,	7,	7,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	9,	9,	9,	9,	9,	9,	9,	9,	9,	9,	10,	10,	10,	10,	10,	10,	10,	10,	10,	10,	10,	10,	10,	10,	11,	11,	11,	11,	11,	11,	11,	11,	11,	11,	11,	11,	11,	11,	11,	11,	12,	12,	12,	12,	12,	12,	12,	12,	12,	12,	12,	12,	12,	12,	12,	12,	12,	12,	12,	12,	12,	12,	12,	12,	13,	13,	13,	13,	13,	13,	13,	13,	13,	13,	13,	13,	13,	13,	13,	13,	14,	14,	14,	14,	14,	14,	14,	14,	14,	14,	14,	14,	14,	14,	14,	14,	15,	15,	15,	15,	15,	15,	15,	15,	15,	15,	16,	16,	16,	16,	16,	16,	16,	16,	16,	16,	16,	16,	16,	16,	16,	16,	16,	16,	16,	16,	17,	17,	17,	17,	17,	17,	17,	17,	17,	17,	17,	17,	17,	17,	17,	17,	17,	17,	17,	17,	18,	18,	18,	18,	18,	18,	18,	18,	18,	18,	18,	18,	18,	18,	18,	18,	18,	18,	18,	18,	18,	18,	18,	18,	18,	18,	19,	19,	19,	19,	19,	19,	19,	19,	19,	19,	19,	19,	19,	19,	19,	19,	19,	19,	19,	19,	19,	19,	19,	19,	19,	19,	19,	19,	19,	19,	19,	19,	20,	20,	20,	20,	20,	20,	20,	20,	20,	20,	20,	20,	20,	20,	20,	20,	20,	20,	21,	21,	21,	21,	21,	21,	21,	21,	21,	21,	21,	21,	21,	21,	21,	21,	21,	21,	22,	22,	22,	22,	22,	22,	22,	22,	22,	22,	22,	22,	22,	22,	22,	22,	22,	22,	22,	22,	22,	22,	22,	22,	22,	22)



#order by network
connectome_whole_360_nodes_network_ordered <- connectome_whole_360_nodes[nodes_network_ordered,nodes_network_ordered]

#rename and cut down number of labels
groups <- c("Primary Visual","", "Early Visual",rep("",5), "Dorsal Stream Visual",rep("",11),"Ventral Stream Visual",rep("",13),"MT+Complex & Neighboring Visual Areas",rep("",17),"Somatosensory & Motor",rep("",9),"Paracentral Lobular & Mid Cingulate",rep("",17),"Premotor",rep("",13),"Posterior Opercular",rep("",9),"Early Auditory",rep("",13),"Auditory Association",rep("",15),"Insular & Frontal Opercular",rep("",23),"Medial Temporal",rep("",15),"Lateral Temporal",rep("",15),"Temporo-Parieto Occipital Junction",rep("",9),"Superior Parietal",rep("",19),"Inferior Parietal",rep("",19),"Posterior Cingulate",rep("",25),"Anterior Cingulate & Medial Prefrontal",rep("",31),"Orbital & Polar Frontal",rep("",17),"Inferior Frontal",rep("",17),"Dorsolateral Prefrontal",rep("",25)) 
rows <- c("Primary Visual","", "Early Visual",rep("",5), "Dorsal Stream Visual",rep("",11),"Ventral Stream Visual",rep("",13),"MT+Complex & Neighboring Visual Areas",rep("",17),"Somatosensory & Motor",rep("",9),"Paracentral Lobular & Mid Cingulate",rep("",17),"Premotor",rep("",13),"Posterior Opercular",rep("",9),"Early Auditory",rep("",13),"Auditory Association",rep("",15),"Insular & Frontal Opercular",rep("",23),"Medial Temporal",rep("",15),"Lateral Temporal",rep("",15),"Temporo-Parieto Occipital Junction",rep("",9),"Superior Parietal",rep("",19),"Inferior Parietal",rep("",19),"Posterior Cingulate",rep("",25),"Anterior Cingulate & Medial Prefrontal",rep("",31),"Orbital & Polar Frontal",rep("",17),"Inferior Frontal",rep("",17),"Dorsolateral Prefrontal",rep("",25)) 

colnames(connectome_whole_360_nodes_network_ordered) <- groups
rownames(connectome_whole_360_nodes_network_ordered) <- rows

#view the connectome on a correlation matrix plot
corrplot(data.matrix(connectome_whole_360_nodes_network_ordered),is.corr=FALSE,tl.srt = 45,tl.cex=0.75,method = "color",col=viridis(200), tl.col = "black")


#if you want to display only sig. values (>.95), delete values which are lower than .95 and replace them with zeros (0) in the matrix.
connectome_whole_sig_360_nodes_network_ordered<-replace(data.matrix(connectome_whole_360_nodes_network_ordered), data.matrix(connectome_whole_360_nodes_network_ordered)<.95, 0) 
colnames(connectome_whole_sig_360_nodes_network_ordered) <- groups

#display sig. nodes only (>.95)
corrplot(connectome_whole_sig_360_nodes_network_ordered,tl.srt = 45,tl.cex=0.75,method = "color",col=viridis(200),is.corr=FALSE,cl.pos="n",tl.col = "black")

#put squares around the networks - can't get it to work :(
#corrplot(connectome_whole_sig_360_nodes_network_ordered,tl.srt = 45,tl.cex=0.75,method = "color",col=viridis(200),is.corr=FALSE,cl.pos="n",tl.col = "black")%>% corrRect(c(1, 100, 150), col="white")




