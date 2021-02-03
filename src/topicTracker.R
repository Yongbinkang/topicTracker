buildTES <- function(x, topic_df, best=c("min","max"), min_tes = -1, ...){
    "
    This function identifies the ancestors of each topic considering its appeared time information.
    Input:
        - x: an evolution strength matrix where a higher value between a past topic and a newer topic indicates that
             the past topic can influence stronger strength on the generation of the new topic. 
        - topic_df: a topic data frame that has all information about topics
        - best: this indicates whether we should choose the max or min when choosing ancestors. But this information
                will be ignore if min_tes is between 0 and 1.
        - min_tes: -1 means that as ancestors we choose the 'best' (latest) one. If this value is in [0,1], 
                    only include the ancestors with evolutionary strengths are >= min_tes
    Output:
        - the tracked topic information consisting of topics and their evolutionary relations over the timeline observed in the topic_df
    "
    ##-----------------------------------------------------------------------------------------------
    ## Check parameters
    best <- match.arg(best)
    if(best=="min"){
        best <- min
        which.best <- which.min
    } else {
        best <- max
        which.best <- which.max
    }
    
    x_names = topic_df$id
    x.dates = topic_df$date
    if(length(x_names) != length(x.dates)){
        stop("inconsistent length for x.dates")
    }

    if(is.character(x.dates)){
        msg <- paste("x.dates is a character vector; " ,
                     "please convert it as dates using 'as.POSIXct'" ,
                     "\n(making sure dates are given as 'YYYY/MM/DD' or 'YYYY-MM-DD').", sep="")
        stop(msg)
    }

    x <- as.matrix(x)
    N <- length(x_names)
    id <- 1:N
    x.dates <- as.POSIXct(round.POSIXt(x.dates,units="days")) # round dates to the day
    
    ## rename dimensions using id
    colnames(x) <- rownames(x) <- id

    if(length(x_names) != nrow(x)){
        stop("inconsistent dimension for x")
    }
    #------------------------------------------------------------------------------------------------

    ##-----------------------------------------------------------------------------------------------    
    ## test equality in floats
    test.equal <- function(val,vec){
        return(abs(val-vec) < 1e-12)
    }
    ##-----------------------------------------------------------------------------------------------

    ##-----------------------------------------------------------------------------------------------
    ## return the names of optimal value(s) in a named vector
    which.is.best <- function(vec){
        res <- names(vec)[test.equal(best(vec), vec)]
        return(res)
    }
    ##-----------------------------------------------------------------------------------------------

    ##-----------------------------------------------------------------------------------------------
    ## select appropriate ancestors among different possible ancestors
    selAmongAncestors <- function(idx, ances, min_tes){
        ances_set = c()
        if (min_tes == -1) {
            # take the latest ancestor
            weight = which.max(x.dates[ances])
            ances <- ances[weight] 
            ances_set[1] = ances
        } else {
            # take the multiple ancestors
            for (i in 1:length(ances)) {
                weight = x[ances[i],idx]
                if (weight >= min_tes) {
                    # If the TES (topic evolution strength) is above min_tes, we add ances[i]
                    ances_set <- append(ances_set, ances[i], after=length(ances_set))
                }
            }
            if (length(ances_set)==0) {
                # if there is no ancestor, by default, we add the ancestor whose weight is the max.
                # This is necessary only for visualisation as we need to make a single tree.
                weight = which.max(x.dates[ances])
                ances <- ances[weight] 
                ances_set <- append(ances_set, ances, after=length(ances_set))
            }
        }
        return(ances_set)
    }
    ##-----------------------------------------------------------------------------------------------

    ##-----------------------------------------------------------------------------------------------
    ## findAncestors
    findAncestors <- function(ids){ 
        entries = hash()
        
        # returns the index of each sequence's ancestor
        res_list = list()
        k = 1
        for (idx in ids) {
            cand_list <- which(x.dates < x.dates[idx])
            if(length(cand_list)==0) {
                # The first generation topics whose ancestor is root
                res_list[[k]] = list(child=idx, ances=NA, weight=0)
                k=k+1
            } else if(length(cand_list)==1) {
                # If there is only one ancestor, we connect the idx to it to make a single tree
                weight = x[cand_list, idx]
                res_list[[k]] = list(child=idx, ances=cand_list, weight=x[cand_list, idx])
                k=k+1
            } else {
                # If there are multiple ancestors, we find best candidates whose TES are >= min_tes.
                anc_list <- selAmongAncestors(idx, cand_list, min_tes) 

                #--------------------------------------------------------------------------
                # If we connect the idx to all ancestors, use the following code block
                "
                for (a in anc_list) {
                    res_list[[k]] = list(child=idx, ances=a, weight=x[a, idx])
                    k = k+1
                }
                "
                #--------------------------------------------------------------------------

                #--------------------------------------------------------------------------
                # If we see mutiple ancestors in the same evolution path, choose the highest one
                removed = c()
                for (a in rev(anc_list)) {
                    if (a %in% removed) next
                    
                    # check if there are super ancestors (ancestors of a (the current ancestor))
                    super_anc_list <- which(x.dates < x.dates[a]) 
                    
                    # remove ancestors whose ids are NA
                    super_anc_list<-super_anc_list[!is.na(super_anc_list)]                                        
                    for (s in super_anc_list) {
                        # If s exists in the current anc_list
                        if (s %in% anc_list) {
                            # Need to compare the weights and choose the one ancestor with the highest weight
                            if (x[a, idx] < x[s, idx]) {
                                key = paste0(idx, s)
                                if (!has.key(key, entries)) {
                                    res_list[[k]] = list(child=idx, ances=s, weight=x[s, idx])
                                    k = k+1    
                                    entries[key] = key
                                }
                            } else {
                                key = paste0(idx, a)
                                if (!has.key(key, entries)) {
                                    res_list[[k]] = list(child=idx, ances=a, weight=x[a, idx])
                                    k = k+1    
                                    entries[key] = key           
                                }
                            }
                            removed <- union(removed, s)
                        } else {
                            key = paste0(idx, a)
                            if (!has.key(key, entries)) {
                                # Add the ancestor for the node - 'a'
                                res_list[[k]] = list(child=idx, ances=a, weight=x[a, idx])
                                k = k+1
                                entries[key] = key                            
                            }
                        }
                    }
                    if (length(super_anc_list)==0) {
                        if (a %in% removed) break
                        key = paste0(idx, a)
                        if (!has.key(key, entries)) {
                            res_list[[k]] = list(child=idx, ances=a, weight=x[a, idx])
                            k = k+1
                            entries[key] = key          
                        }
                    }
                 }
                #--------------------------------------------------------------------------
            }
        }
        res = do.call(cbind, res_list)
        return(res)
    }
    ##-------------------------------------------------------------------------


    ## BUILD THE OUTPUT ##
    res <- findAncestors(id)
    res <- data.frame(child=unlist(res[1,]), ances=unlist(res[2,]), weight=unlist(res[3,]))
    ances.date <- x.dates[res[,2]]
    child.name <- x_names[res[,1]]
    ances.name <- x_names[res[,2]]
    label <- topic_df$label[res[,1]]
    child.word <- topic_df$words[res[,1]]
    
    res <- data.frame(child.name=unlist(child.name), ances.name=unlist(ances.name), res, 
                      child.date=x.dates[res[,1]], ances.date, label, child.word)
    class(res) <- c("topicTrack","data.frame")
    
    return(res)
} 

## -------------------------------------------------------------
## These functions are utilise functions for visualiseTES()

t_col <- function(color, percent = 10, name = NULL) {
    "
    This function makes the transparent color of the input color. This function will be used to color 
    the legend.
    - color: input color name
    - percent: transparent ratio (% transparency
    - name = an optional name for the color
    "

    ## Get RGB values for named color
    rgb.val <- col2rgb(color)
    ## Make new color using input color as base and alpha set by transparency
    t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
           max = 255,
           alpha = (100-percent)*255/100,
           names = name)
    ## Save the color
    invisible(t.col)
}

update_mark_group <- function (topics, x_names, evol_color_group, mark_colors, 
                               legend_colors, legend_color_index, mark_expand, k, mark_size){
    "
    This function updates the colors of the nodes based on their evolutionary status
    "
    indices = which(x_names %in% topics)
    for (index in indices) {
        evol_color_group[[k]] = index
        mark_colors[k] = legend_colors[legend_color_index]
        k = k+1
    }      
    for (i in length(mark_expand):k) {
        mark_expand[i] = mark_size
    }
    
    return (list(evol_color_group=evol_color_group, 
                 mark_colors=mark_colors, 
                 mark_expand=mark_expand, 
                 k=k))
}
## -------------------------------------------------------------

visualiseTES <- function(x, topic_df, root='R', min_reborn=2, min_dead=2, min_tes=0.2, out_image_path, ...){

    "
    The function taking care of visualisation of topic trajectory information in TET (topic evolution tree).
    Input:
    - x: the topictrack data frame that is the result of the 'topicTrack' function.
    - root: label of the root node of the TET
    - min_reborn: the minimum distance between infector and infectee nodes
    - min_dead: the minimum sleeping years
    - min_tes: the minimum of the topic evolution strength (test) to be visualised
    "

    #---------------------------------------------------
    # Read the number of nodes to be visualised except root
    nr_of_hosts <- length(unique(x$child))
    rootnodetime = min(x$child.date)
    tailnodetime = max(x$child.date)
    rootnodetime_year <- as.numeric(format(rootnodetime,'%Y'))
    tailnodetime_year <- as.numeric(format(tailnodetime,'%Y'))
    labs <- pretty(range(par("xaxp")[1:2] + rootnodetime, par("xaxp")[1:2] + tailnodetime))
    pos = labs - (rootnodetime)
    #---------------------------------------------------
    
    #---------------------------------------------------
    # Create a topic evolution tree (TET) using layout_as_tree
    x$ances[is.na(x$ances)] <- 0
    # If the node with weight is less than min_tes, we connect the node to the root.
    x$weight[x$weight<=min_tes] <- 0
    x$ances[(x$ances !=0) & (x$weight == 0)] <- 0
    
    topic_track_df <- data.frame(x=x$ances, y=x$child, weight=x$weight)
    G <- graph.data.frame(d=topic_track_df, directed = TRUE)
    hosts = unique(as.vector(x$child.name))
    hosts_with_root = c(root, hosts)
    x_names = (hosts_with_root)[1+unique(c(topic_track_df[,1], 1:nr_of_hosts))]
    x_names_without_root = x_names[2:length(x_names)]
    layout_G = layout_as_tree(G)[,c(2,1)]*rep(c(-1,1),each=1+nr_of_hosts)
    #-----------------------------------------

    #---------------------------------------------------
    # colours for the arrows
    edgeweight <- x$weight
    edgeweight[is.na(edgeweight)] <- 0
    ewreds <- (edgeweight > 0.4) -.5*(edgeweight > 0.8)
    ewgreens <- .5*(edgeweight > 0.2) -.5*(edgeweight > 0.6)
    ewblues <- (edgeweight <= 0.2) + .5*(edgeweight > 0.8)
    #---------------------------------------------------
    
    #---------------------------------------------------
    # Get the node weights
    topic_weights = c()
    topic_weights[1] = 0
    i = 2
    for (name in x_names_without_root) {
        topic_weights[i] = topic_df[which(topic_df$id==name)[1],]$weight
        i =i+1
    }
    
    # Sort the topics according to their weights in y-axis
    y_sort = rank(topic_weights)
    for (i in 1:length(topic_weights)) {
        layout_G[,2][i] = topic_weights[i]
    }
    #---------------------------------------------------
    
    #---------------------------------------------------    
    # Set names of the nodes in the graph
    i = 1
    # Get labels of topics corresponding to x_names
    x_labels = c()
    for (name in x_names_without_root) {
        id = as.character(x[which(x$child.name==name),]$label[1])
        label = as.character(x[which(x$child.name==name),]$label[1])
        x_labels[i] = name
        i = i+1
    }
    x_labels = c(root, x_labels)
    
    # Read positions of the nodes in the graph
    x_vec = layout_G[,1]
    y_vec = layout_G[,2]
    #---------------------------------------------------    
    
    #---------------------------------------------------
    # Adjust the (x,y) coords of all nodes based on their dates (years)
    
    # First, X-coord fitting: the nodes with the same year must have the same x-coord.
    # Read the year information of the nodes.
    x_min = min(x_vec)
    x_max = max(x_vec)
    year_vectors = c()
    year_vectors[1] = rootnodetime_year-1
    i = 2
    for (name in x_names_without_root) {
        year_vectors[i] = as.numeric(format(x[which(x$child.name==name)[1],]$child.date,'%Y'))
        i =i+1
    }
    
    # Adjust x-coords of all nodes
    y_dict = hash()
    y_index_dict = hash()
    
    x_vec_adjusted = x_vec
    year_min = min(year_vectors)
    year_max = max(year_vectors)
    for (i in 2:length(x_vec)) {
        year = year_vectors[i]
        key = toString(year)

        # change x-coords in the way that all nodes with the same year have the same x-coord
        x_current = x_vec[i]
        x_truth = (year-year_min) + x_min
        x_vec_adjusted[i] = x_truth
        layout_G[,1][i] = x_truth
        
        # store y-coords for resolving duplicated position in the next step
        if (is.null(y_dict[[key]])) {
            y_dict[[key]] = round(y_vec[i], 4)
            y_index_dict[[key]] = i
        } else {
            y_dict[[key]] = c(y_dict[[key]], round(y_vec[i], 4))
            y_index_dict[[key]] = c(y_index_dict[[key]], i)
        }        
    }
    # ----------------------------------------------------------------------------
    
    # ----------------------------------------------------------------------------
    # Group the nodes according to their dates. Let's first make the group of the colors to be used.
    k = 1
    color_groups = list()
    unique_xcoords = unique(x_vec_adjusted)
    for (i in 2:length(unique_xcoords)) {
        group = c(which(unique_xcoords[i]==x_vec_adjusted))
        color_groups[[k]] = group
        k = k+1
    }
    valcol <- hcl.colors(100, "heat2", rev = TRUE)
    # ----------------------------------------------------------------------------    

    # ----------------------------------------------------------------------------
    # Group the nodes according to their types: These information will be used to grouping in colors
    # 1) newly born: 
    # 2) fused node: if a node is evolved from more than 2 previous nodes.
    # 3) split a node: if a node is split to next year's nodes.
    # 4) dead node: if a node has no infectees before a user-specified period - min_dead
    # 5) born again: if a node emerges after a user-specified period - min_reborn
    # 6) flourishing: by default all nodes are flourishing except the above cases.
    
    k = 1
    legend_color_index = 1
    evol_color_group = list()
    mark_colors = c()
    legend_colors = c(t_col("gold1", percent = 10),
                      t_col("orange", percent = 10),
                      t_col("cyan", percent = 10),
                      t_col("gray", percent = 10),
                      t_col("yellowgreen", percent = 10)
                      )
    #legend_colors = rainbow(5, alpha=.1)
    mark_expand=c()
    
    ### Get newly born nodes
    # if you want to change the node size according its evolution state, please change the value of mark_size below

    children_of_root = (x[which(topic_track_df$x==0),])[,1]
    children_of_weak_nodes = (x[which(x$weight==0),])[,1]
    topics = unique((union(children_of_root, children_of_weak_nodes)))
    res = update_mark_group (topics, x_names, evol_color_group, mark_colors, 
                               legend_colors, legend_color_index, mark_expand, k, mark_size=20)
    legend_color_index = legend_color_index+1
    
    ### Get fused nodes
    df <- data.frame(x=x$ances, y=x$child)
    # Get node indices of duplicated infectees
    node_indices = df[duplicated(df$y),2]
    # Get topics names corresponding to the node_indices
    topics = unique(x[which(x$child %in% node_indices),][,1])
    res = update_mark_group (topics, x_names, res$evol_color_group, res$mark_colors, 
                               legend_colors, legend_color_index, res$mark_expand, res$k, mark_size=15)
    legend_color_index = legend_color_index+1
    
    ### Get split nodes
    # Get node indices of duplicated infectors
    node_indices = df[duplicated(df$x),1]
    # Get topics names corresponding to the node_indices
    topics = unique(x[which(x$child %in% node_indices),][,1])
    res = update_mark_group (topics, x_names, res$evol_color_group, res$mark_colors, 
                               legend_colors, legend_color_index, res$mark_expand, res$k, mark_size=10)
    legend_color_index = legend_color_index+1
    
    ### Get dead nodes
    # Get node indices of all infectees
    node_indices = unique(topic_track_df[2])[,1]
    # Get topics names corresponding to the node_indices
    topics = x[which(x$child %in% node_indices),]
    # Find the nodes that do not have children (i.e. outgoing edges)
    child_ids_without_outgoing = unique(setdiff(topics$child,topics$ances))
    child_ids_without_outgoing_indices = which(topics$child %in% child_ids_without_outgoing)
    topics = topics[child_ids_without_outgoing_indices,]
    # step1: Get the topics where their born dates are different from the latest year
    max_date = as.numeric(format(max(x$child.date),'%Y')) 
    child_dates = as.numeric(format(topics$child.date,'%Y'))    
    topics = unique(topics[which(max_date-child_dates > min_dead),])[,1]
    res = update_mark_group(topics, x_names, res$evol_color_group, res$mark_colors, 
                               legend_colors, legend_color_index, res$mark_expand, res$k, mark_size=5)
    legend_color_index = legend_color_index+1
    
    ### Get reborn nodes
    node_indices = unique(df[2])[,1]
    # Get topics names corresponding to the node_indices
    topics = x[which(x$child %in% node_indices),]
    # Get the topics whose date difference with their infectors is more than 'min_reborn' years
    latest = unique(subset(topics, ances.date==ave(ances.date, child, FUN=max)))
    child_dates = as.numeric(format(latest$child.date,'%Y')) 
    ances_dates = as.numeric(format(latest$ances.date,'%Y'))    
    topics = latest[which(child_dates-ances_dates > min_reborn),][,1]    
    res = update_mark_group(topics, x_names, res$evol_color_group, res$mark_colors, 
                               legend_colors, legend_color_index, res$mark_expand, res$k, mark_size=3)
    legend_color_index = legend_color_index+1
    
    evol_color_group = res$evol_color_group
    mark_colors = res$mark_colors
    mark_expand = res$mark_expand
    # ----------------------------------------------------------------------------
    
    # ----------------------------------------------------------------------------
    # Plot the graph
    
    y_min = min(y_vec)
    y_max = max(y_vec)
    layout_G = norm_coords(layout_G, ymin=-1, ymax=1, xmin=-1, xmax=1)
    
    vertex_size = c(0,rep(3, nr_of_hosts))[1+unique(c(topic_track_df[,1],1:nr_of_hosts))]
    plot(G, layout=layout_G,
         #vertex.size = vertex_size, 
         vertex.size = 10, 
         edge.arrow.size = .3,
         vertex.label=NA,
         vertex.label.cex = .7, 
         vertex.shape="circle",
         #vertex.label.degree = -pi/2,
         #vertex.label = ifelse(x_names==root, "", x_names),
         #vertex.label = ifelse(x_names==root, "", x_labels),         
         vertex.label.color = "black",
         vertex.color = NA, 
         #vertex.frame.color = "black",
         vertex.frame.color = ifelse(x_names==root, "white", "black"),
         vertex.frame.width = 5, 
         #edge.color = rgb(ewreds,ewgreens,ewblues), 
         edge.color = ifelse(ewblues==1, "NA", rgb(ewreds,ewgreens,ewblues)), 
#          edge.curved=.3,
         edge.width = 0.5, axes = F, edge.curved=F,
         vertex.label.dist=0, 
         #mark.groups=color_groups, mark.border = NA, mark.col=valcol
         mark.groups=evol_color_group, mark.border = NA, mark.col=mark_colors, mark.expand=mark_expand,
         alpha = .1,
         xlim=c(-1,1),
         ylim=c(-1,1),
         rescale=F,
         xlab="Timeline (Years)",
         ylab="Topic weights",
         asp=1
        )    

    # To avoid the overlapped labels of the nodes
    thigmophobe.labels(layout_G[,1], layout_G[,2], labels=ifelse(x_names==root, "", x_labels), 
                       #font=0.5,
                       cex=0.75, offset=0.8)
    
    legend("bottomright", 
           title="Evolution strength",
           legend=c("[20%,40%)", "[40%,60%)", "[60%,80%)", "[80%,100%)"),
           col=c("darkgreen", "orange", "red", "darkviolet"), 
           cex=0.8,
           lty=1, 
           text.font=0.5, bg='white',
           xpd=TRUE,
           horiz=F,
           inset=c(-.3,0),
           bty="n",
          ) 
    
    legend("topright", 
           title="Evolution states",
           #inset=.02, 
           #inset=c(0, -.15),
           legend=c("born", "fused", "split", "dead", "reborn"),
           col=legend_colors, 
           cex=0.8,
           #box.lty=1, box.lwd=1, box.col="black", 
           text.font=0.5, 
           bg='white',
           pch=c(16,16,16,16,16),
           pt.cex = 2,
           xpd=TRUE,
           horiz=F,
           inset=c(-.3,0),
           bty="n"
          )
    # ----------------------------------------------------------------------------
    

    # ----------------------------------------------------------------------------
    # Draw the ticks
    
    # Define the ticks of x-axis (timeline)
    no_ticks = tailnodetime_year-rootnodetime_year+1
    no_ticks = no_ticks+1
    ticks_by = round(2/(no_ticks), digits=2)
    
    pos = seq(from=-1, to=1, length.out=no_ticks)

    print_labs = c()
    print_labs[1] = rootnodetime_year-1
    for (i in 2:(no_ticks)) { 
        print_labs[i]=1+print_labs[i-1]
    }
    axis(side=1, at=pos, labels=print_labs, cex.axis=0.8)

    y_no_ticks = 5
    ytick <- seq(-1, 1, 0.5)
    tmp <- seq(y_min, y_max, (y_max-y_min)/4)
    y_labs = c()
    for (i in 1:(y_no_ticks)) { 
        y_labs[i] = round(tmp[i], 2)
    }
    
    axis(side=2, at=ytick, labels=y_labs, cex.axis=0.8)
    
    # ----------------------------------------------------------------------------
    result <- list("G"=G, "layout"=layout_G, "node_ids"=x_names, "node_labels"=x_labels, "node_size")
    return(result)
}

