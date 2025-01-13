scalers = {
    'c': "<operator spec='ScaleOperator' id='cScaler' parameter='@c' scaleFactor='0.80' weight='1'/>",
    'K': "<operator spec='ScaleOperator' id='KScaler' parameter='@K' scaleFactor='0.80' weight='1'/>",
    'theta': "<operator spec='ScaleOperator' id='thetaScaler' parameter='@theta' scaleFactor='0.80' weight='1'/>",
    'alpha': "<operator spec='ScaleOperator' id='alphaScaler' parameter='@alpha' scaleFactor='0.80' weight='1'/>",
    'mu': "<operator spec='ScaleOperator' id='muScaler' parameter='@mutationRate' scaleFactor='0.80' weight='1'/>"
}	



c_priors = {
    1.0: "<distr spec='LogNormalDistributionModel' M='0' S='0.5'/>",
    10.0: "<distr spec='LogNormalDistributionModel' M='2.5' S='0.5'/>"
}

K_priors = {
    0.1: "<distr spec='LogNormalDistributionModel' M='-2' S='0.5'/>",
    1.0: "<distr spec='LogNormalDistributionModel' M='0' S='0.5'/>",
    10.0: "<distr spec='LogNormalDistributionModel' M='2.5' S='0.5'/>"
}

theta_priors = {
    1.0: "<distr spec='LogNormalDistributionModel' M='0.0' S='0.5'/>"
}

alpha_priors = {
    0.1: "<distr spec='Beta' alpha='1' beta='9'/>",
    0.5: "<distr spec='Beta' alpha='10' beta='10'/>",
    0.99: "<distr spec='Beta' alpha='9' beta='1'/>"
}

mu_priors = {
    0.01: "<distr spec='LogNormalDistributionModel' M='-2.0' S='2.0'/>"
}

kappa_priors = {
    3.0: "<distr spec='LogNormalDistributionModel' M='1.0' S='0.5'/>"
}


inital_tree_constant = {
    False: "<init spec='SeedbankTreeInitializer'>\
			\n\t\t\t<transitionModel idref='transitionModel'/>\
			\n\t\t\t<initial idref='sbTree'/>\
			\n\t\t\t<trait idref='typeTraitSet'/>\
			\n\t\t\t<trait idref='timeTraitSet'/>\
		\n\t\t</init>",
    True: "<init spec='SeedbankTreeFromNewick' \
		\n\t\tvalue='$$$REPLACE_NEWICK$$$'\
		\n\t\t> \
			\n\t\t\t<initial idref='sbTree'/>\
			\n\t\t\t<data idref='alignment'/>\
		\n\t\t</init>",
}


tree_operators = "<operator spec='SeedbankTreeScale' id='treeScale' sbTree='@sbTree' transitionModel='@transitionModel' scaleFactor='0.80' weight='1'/>\
\n\t\t<operator spec='TypedWilsonBalding' sbTree='@sbTree' transitionModel='@transitionModel' weight='10'/>\
\n\t\t<operator spec='TypedSubtreeExchange' sbTree='@sbTree' transitionModel='@transitionModel' weight='10' isNarrow='true'/>\
\n\t\t<operator spec='TypedSubtreeExchange' sbTree='@sbTree' transitionModel='@transitionModel' weight='10' isNarrow='false'/>\
\n\t\t<operator spec='NodeShiftRetype' sbTree='@sbTree' transitionModel='@transitionModel' weight='10' rootOnly='true' rootScaleFactor='0.8'/>\
\n\t\t<operator spec='NodeShiftRetype' sbTree='@sbTree' transitionModel='@transitionModel' weight='10' rootOnly='false'/>"