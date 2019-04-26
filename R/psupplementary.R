#' Dimensionality reduction plots and standard psupertime results for acinar cells
#' (Fig 1C, Fig 1D, Fig 1E, Supp Fig 01, Supp Fig 02, Supp Fig 03)
#'
#' @import data.table
#' @import ggplot2
#' @export
acinar_cells_plots <- function(output_dir='.', ext='png') {
	# check whether can use umap
	if ( !requireNamespace("umap", quietly=TRUE) ) {
		message('umap not installed; not doing clustering')
		return()
	}
	library('umap')

	# get data
	message('loading acinar cell SCE object')
	tag 			= 'aging_acinar'
	label_name 		= 'Donor age\n(years)'
	data(acinar_sce)

	# restrict to highly variable genes
	message('identifying highly variable genes')
	hvg_params 		= list(hvg_cutoff=0.1, bio_cutoff=0.5, span=0.1)
	sel_genes 		= psupertime:::calc_hvg_genes(acinar_sce, hvg_params)
	acinar_hvg 		= acinar_sce[sel_genes, ]

	# calc and plot umap
	message('projecting using UMAP')
	x 				= t(SummarizedExperiment::assay(acinar_hvg, 'logcounts'))
	wellKey_vector 	= SingleCellExperiment::colData(acinar_hvg)$wellKey
	label_vector 	= factor(SingleCellExperiment::colData(acinar_hvg)[['donor_age']])
	proj_umap 		= calc_umap(x, wellKey_vector)
	g 				= plot_dim_reduction(label_vector, proj_umap$umap1, proj_umap$umap2, labels=c('umap1', 'umap2', 'Donor age (years)'))
	plot_file 		= file.path(output_dir, sprintf('fig1C %s_umap_plot.%s', tag, ext))
	ggplot2::ggsave(plot_file, g, height=6, width=7)

	# run psupertime
	message('running psupertime')
	y_age 			= acinar_sce$donor_age
	psuper_obj 		= psupertime(acinar_sce, y_age)

	# plot labels over learned psupertime
	message('plotting labels over psupertime')
	g 				= plot_labels_over_psupertime(psuper_obj, label_name)
	plot_file 		= file.path(output_dir, sprintf('fig1D %s labels over psupertime.%s', tag, ext))
	ggplot2::ggsave(plot_file, g, height=6, width=12)

	# plot identified genes against learned psupertime
	message('plotting identified genes over psupertime')
	g 				= plot_identified_genes_over_psupertime(psuper_obj, label_name, n_to_plot=5)
	plot_file 		= file.path(output_dir, sprintf('fig1E %s labels over psupertime.%s', tag, ext))
	ggplot2::ggsave(plot_file, g, height=2, width=10)
	g 				= plot_identified_genes_over_psupertime(psuper_obj, label_name)
	plot_file 		= file.path(output_dir, sprintf('suppfig02 %s identified genes over psupertime.%s', tag, ext))
	ggplot2::ggsave(plot_file, g, height=8, width=12)

	# plot training diagnostics
	message('plotting psupertime training diagnostics')
	g 				= plot_train_results(psuper_obj)
	plot_file 		= file.path(output_dir, sprintf('suppfig01 %s training results.%s', tag, ext))
	ggplot2::ggsave(plot_file, g, height=6, width=6)

	# plot coefficients for identified genes
	message('plotting identified gene coefficients')
	g 				= plot_identified_gene_coefficients(psuper_obj)
	plot_file 		= file.path(output_dir, sprintf('suppfig03 %s identified genes.%s', tag, ext))
	ggplot2::ggsave(plot_file, g, height=6, width=8)
}

#' Clustering of genes, plus GO term analysis of them
#' (Supp Fig 04, Supp Fig 05)
#'
#' @import data.table
#' @import ggplot2
#' @export
acinar_cells_go_term_analysis <- function(output_dir='.', ext='png') {
	# get data
	tag 			= 'aging_acinar'
	label_name 		= 'Donor age\n(years)'
	data(acinar_sce)

	# run psupertime
	y_age 			= acinar_sce$donor_age
	psuper_obj 		= psupertime(acinar_sce, y_age)

	# calculate go terms
	go_list 		= psupertime_go_analysis(psuper_obj, org_mapping='org.Hs.eg.db')

	# plot gene cluster profiles
	g 				= plot_profiles_of_gene_clusters(go_list, label_name='Donor age\n(years)')
	plot_file 		= file.path(output_dir, sprintf('suppfig04 %s gene cluster profiles.%s', tag, ext))
	ggsave(plot_file, g, height=10, width=8)

	# plot go terms
	g 				= plot_go_results(go_list)
	plot_file 		= file.path(output_dir, sprintf('suppfig05 %s significant GO terms by cluster.%s', tag, ext))
	ggsave(plot_file, g, height=12, width=6)
}

#' Calculates UMAP projection
#'
#' @import data.table
#' @internal
calc_umap <- function(x, wellKey_vector) {
	set.seed(1)
	cat('calculating umap\n')
	umap_obj 	= umap::umap(x)
	proj_umap 	= data.table(
		wellKey 	= wellKey_vector,
		umap_obj$layout
		)

	# scale
	setnames(proj_umap, c('V1', 'V2'), c('umap1', 'umap2'))
	proj_umap[, umap1 := scale_dim(umap1) ]
	proj_umap[, umap2 := scale_dim(umap2) ]

	return(proj_umap)
}

#' Calculates tSNE projection
#'
#' @import data.table
#' @internal
calc_tsne <- function(x, wellKey_vector) {
	set.seed(1)

	# make wide, do PCA on these
	cat('calculating tsne\n')
	tsne_obj  	= Rtsne(x, initial_dims=20, perplexity=30)
	proj_tsne 	= data.table(
		wellKey 	= wellKey_vector,
		tsne_obj$Y
		)

	# scale
	setnames(proj_tsne, c('V1', 'V2'), c('tsne1', 'tsne2'))
	proj_tsne[, tsne1 := scale_dim(tsne1) ]
	proj_tsne[, tsne2 := scale_dim(tsne2) ]

	return(proj_tsne)
}


#' Scales dimensions for nice plotting
#'
#' @internal
scale_dim <- function(v) {
	v_range 	= range(v)
	v_out 		= (v - v_range[1])*0.9 / (v_range[2] - v_range[1]) + 0.05

	return(v_out)
}

#' Nice plotting of dimensionality reduction
#'
#' @import data.table
#' @import ggplot2
#' @internal
plot_dim_reduction <- function(y_vector, dim1, dim2, labels) {
	# plot pca
	plot_dt 	= data.table(
		y_var 		= y_vector
		,dim1 		= dim1
		,dim2 		= dim2
		)
	col_vals 	= psupertime:::make_col_vals(plot_dt$y_var)
	g = ggplot(plot_dt) +
		aes( x=dim1, y=dim2, colour=y_var ) +
		geom_point() +
		scale_colour_manual( values=col_vals ) +
		# scale_x_continuous( breaks=pretty_breaks(), limits=c(0,1) ) +
		# scale_y_continuous( breaks=pretty_breaks(), limits=c(0,1) ) +
		labs(
			x 		= labels[[1]],
			y 		= labels[[2]],
			colour 	= labels[[3]]
			) +
		theme_bw() +
		theme(
			axis.text 	= element_blank()
			)

	return(g)
}

#' Clustering and psupertime analysis of colon cells
#' (Fig 1F, Supp Fig 12, Supp Fig 13, Supp Fig 14, Supp Fig 15)
#'
#' @import data.table
#' @import ggplot2
#' @export
unsupervised_clustering_of_colon_cells <- function(output_dir='.', ext='png') {
	# is Seurat installed?
	if ( !requireNamespace("Seurat", quietly=TRUE) ) {
		message('Seurat not installed; not doing clustering')
		return()
	}
	if ( !requireNamespace("umap", quietly=TRUE) ) {
		message('umap not installed; not doing clustering')
		return()
	}
	if ( !requireNamespace("cowplot", quietly=TRUE) ) {
		message('cowplot not installed; not doing clustering')
		return()
	}

	# get data
	data(colon_sce)

	# calculate HVGs
	message('calculating highly variable genes')
	hvg_params 		= list(hvg_cutoff=0.5, bio_cutoff=0, span=0.01)
	sel_genes 		= psupertime:::calc_hvg_genes(colon_sce, hvg_params)
	sce_hvg 		= colon_sce[sel_genes, ]

	# do umap, get layout from it
	message('doing UMAP projection')
	set.seed(1)
	umap_obj 		= umap::umap(t(SummarizedExperiment::assay(sce_hvg, 'logcounts')))
	layout_dt 		= data.table(umap_obj$layout)
	setnames(layout_dt, names(layout_dt), c('UMAP1', 'UMAP2'))

	# do clustering with Seurat
	message('unsupervised clustering with Seurat')
	clusters 		= do_seurat_clustering(sce_hvg)

	# add clusters, put in sensible order
	layout_dt[, raw_cluster := clusters ]
	cluster_order 	= data.table(
		raw_cluster = factor(as.character(c(0, 1, 2, 3, 4, 5, 6, 7, 8))),
		cluster 	= factor(as.character(as.integer(c(8, 4, 3, 6, 9, 1, 2, 5, 7))))
		)
	# cluster_order 	= layout_dt[, list(cluster_med = median(UMAP1)), by=raw_cluster]
	# cluster_order[, cluster := factor(as.character(rank(cluster_med))) ]
	layout_dt 		= cluster_order[layout_dt, on='raw_cluster']

	# plot unlabelled
	message('plotting clusterings')
	g = ggplot(layout_dt) +
		aes( x=UMAP1, y=UMAP2) +
		geom_point(size=3, fill='grey', colour='black', shape=21 ) +
		theme_light() +
		theme(
			axis.text 	= element_blank()
			)
	umap_file 		= file.path(output_dir, sprintf('fig1F colon unlabelled hvg umap.%s', ext))
	ggsave(umap_file, g, height=6, width=8)

	# plot unsupervised clusters
	g = ggplot(layout_dt) +
		aes( x=UMAP1, y=UMAP2, colour=cluster) +
		geom_point() +
		scale_colour_brewer( palette='Set1' ) +
		theme_light() +
		theme(
			axis.text 	= element_blank()
			)
	umap_file 		= file.path(output_dir, sprintf('fig1F colon unsupervised hvg umap.%s', ext))
	ggsave(umap_file, g, height=6, width=8)

	# define things we want to do
	order_list 		= list(
		c('1', '4', '6', '8'),
		c('1', '2', '3', '5', '9')
		)
	palette 		= 'BrBG'
	label_name 		= 'Selected\nSeurat\nclusters'
	g_list 			= list()
	ii  			= 1
	base_size 		= 12

	# do psupertime on different cluster sequences
	message('analysing selected cluster orderings')
	for (this_order in order_list) {
		# define label
		tag 			= paste0(this_order, collapse='')
		message('analysing sequence ', tag)
	
		# do psupertime
		y_test_all 		= factor(layout_dt$cluster, levels=this_order)
		keep_idx 		= !is.na(y_test_all)
		sce_test 		= sce_hvg[, keep_idx]
		y_test 			= y_test_all[ keep_idx ]
	
		# do psupertime
		message('running psupertime on this sequence')
		psuper_obj 		= psupertime(sce_test, y_test, sel_genes='all')

		# plot genes identified
		message('plotting genes over identified psupertime')
		g 				= plot_identified_genes_over_psupertime(psuper_obj, label_name='Ordered\nSeurat\nclusters', palette=palette, n_to_plot=20)
		plot_file 		= file.path(output_dir, sprintf('suppfig17_18 %s identified genes over psupertime.%s', tag, ext))
		ggplot2::ggsave(plot_file, g, height=8, width=12)

		# plot selected cluster ordering
		message('plotting selected sequence over UMAP')
		g 				= plot_selected_cluster_ordering(layout_dt, this_order)
		plot_file 		= file.path(output_dir, sprintf('fig1F unsupervised cluster ordering %s.%s', tag, ext))
		ggsave(plot_file, g, height=4, width=5)

		# do go stuff
		message('GO enrichment analysis of gene clusters')
		go_list 		= psupertime_go_analysis(psuper_obj, org_mapping='org.Mm.eg.db')
			
		# assemble supplementary plots
		message('plotting enrichment results')
		g_list[[ii]] 	= plot_profiles_of_gene_clusters(go_list, label_name=label_name, palette=palette)
		ii 				= ii + 1
		g_list[[ii]] 	= plot_go_results(go_list) + theme( axis.text=element_text(size=6))
		ii 				= ii + 1
	}

	# plot supplementary
	message('larger plot of enrichment results')
	g 			= cowplot::plot_grid(plotlist=g_list, labels=LETTERS[1:4], ncol=2, align='h', axis='b')
	plot_file 	= file.path(output_dir, sprintf('suppfig16 colon go terms.%s', ext))
	cowplot::save_plot(plot_file, g, ncol=2, base_height=12, base_aspect_ratio=1/3)
}

# Unsupervised clustering
#'
#' @internal
do_seurat_clustering <- function(sce_hvg) {
	# set up seurat
	seurat_obj 	= Seurat::CreateSeuratObject(
		SummarizedExperiment::assay(sce_hvg, 'counts'), project="temp", 
		min.cells=0, min.genes=0
		)
	seurat_obj@scale.data 	= SummarizedExperiment::assay(sce_hvg, 'logcounts')

	# do PCA and clustering
	seurat_obj 	= Seurat::RunPCA(object=seurat_obj, pc.genes=rownames(sce_hvg), pcs.compute=20, do.print=FALSE)
	seurat_obj 	= Seurat::FindClusters(
		object=seurat_obj, 
		reduction.type="pca", dims.use=1:6, 
		resolution=0.6, print.output=0, save.SNN=TRUE
		)
	clusters 	= seurat_obj@ident

	return(clusters)
}

#' Plot this clustering nicely
#' @import data.table
#' @import ggplot2
#' @internal
plot_selected_cluster_ordering <- function(layout_dt, this_order) {
	palette 		= 'BrBG'

	# get clusters to exclude
	grey_clusts 	= setdiff(layout_dt$cluster, this_order)
	n_ordered 		= length(this_order)
	n_grey 			= length(grey_clusts)
	ordered_vals 	= rev(RColorBrewer::brewer.pal(n_ordered, palette))
	grey_vals 		= RColorBrewer::brewer.pal(n_grey+1, 'Greys')[-1]
	pal_vals 		= c(ordered_vals, grey_vals)
	sort_idx 		= order(c(this_order, grey_clusts))

	# plot results
	g = ggplot(layout_dt) +
		aes( x=UMAP1, y=UMAP2, fill=cluster) +
		geom_point(size=3, colour='black', shape=21 ) +
		scale_fill_manual( values=pal_vals[sort_idx] ) +
		theme_light() +
		theme( axis.text=element_blank() ) +
		labs(
			fill 		= 'Cluster'
			)
	return(g)
}
