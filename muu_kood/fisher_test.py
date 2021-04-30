import scipy.stats as stats

#Nullhüpotees - RNA splaissimise ja geeniekspressiooni puhul on ühe geeni ülekatete suhe mitmega sama
#kui p-väärtus on väiksem kui 0,05, siis saame nullhüpoteesi kummutada

oddsratio, pvalue = stats.fisher_exact([[33511, 7807], [19896, 1701]]) #cs = 50
print(pvalue) #kas on statistiliselt oluline (statistiline seos)

