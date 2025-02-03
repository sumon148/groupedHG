# Install The package
devtools::check()
library(groupedHG)
?pmfHG.perfect
pmfHG.perfect(ty = 0, N = 100, barN = 4, Tx = 20, b = 4)

?pmfHG.imperfect.group
pmfHG.imperfect.group(ty = 0, N = 100, barN = 4, Tx = 20, b = 4, delta = 1, lambda =1)
pmfHG.imperfect.group(ty = 0, N = 100, barN = 4, Tx = 20, b = 4, delta = 0.70, lambda =1)
pmfHG.imperfect.group(ty = 0, N = 100, barN = 4, Tx = 20, b = 4, delta = 0.70, lambda =0.80)

?pmfHG.imperfect.item
pmfHG.imperfect.item(ty = 0, N = 100, barN = 4, Tx = 20, b = 4, delta = 1, lambda =1)
pmfHG.imperfect.item(ty = 0, N = 100, barN = 4, Tx = 20, b = 4, delta = 0.70, lambda =1)
pmfHG.imperfect.item(ty = 0, N = 100, barN = 4, Tx = 20, b = 4, delta = 0.70, lambda =0.80)


?qDLkGroup
?dqDLkGroup
?qdlkItem
?dqdlkItem


?ExpVarBN.imperfect.group

ExpVarBN.imperfect(N = 1000, Tx = 10, barN = 50, b = 5, delta = 0.1, lambda = 0.05, type = "group")
ExpVarBN.imperfect(N = 1000, Tx = 10, barN = 50, b = 5, delta = 0.1, lambda = 0.05, type = "item")
