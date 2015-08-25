#UI

shinyUI( fluidPage(
	titlePanel('Hierarchial Clustering Plot'),

	sidebarPanel(
		#Select Step
		selectInput('step', label = h4('Step'),
			choices = list('Raw Data' = 'raw', 'Normalized' = 'norm', 'Average Rep Probes' = 'exprs'),
			selected = 'norm'),

		#Select Distance Method
		selectInput('dist', label = h4('Distance Method:'),
			choices = list('Euclidean' = "euclidean",'Maximum' = "maximum",'Manhattan' = "manhattan", 
			'Canberra' = "canberra",'Binary' = "binary",'Minkowski' = "minkowski"),
			selected = "euclidean"),
		#Select Clustering Method
		selectInput('clust', label = h4('Clust method:'),
			choices = list('Ward.D' = "ward.D",'Ward.D2' = "ward.D2",'Single' = "single",
			'Complete' = "complete",'Average' = "average", 
			'Mcquitty' = "mcquitty",'Median' = "median",'Centroid' = "centroid"),
			selected = "single"),
		helpText('Average=UPGMA, Mcquitty=WPGMA, Median=WPCMC, Centroid=UPGMC')
	),
	mainPanel(
		#OUTPUTS
		plotOutput('plot')
	)
))