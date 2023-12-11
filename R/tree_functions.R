# Creating a stump for a tree
stump <- function(data,
                  j){

  # Creating the base node
  node <- list()


  node[["node0"]] <- list(
    # Creating the node number
    node_number = 0,
    j = j,
    pred_vars = j,
    inter = NA,
    isRoot = TRUE,
    # Creating a vector with the tranining index
    train_index = 1:nrow(data$x_train),
    test_index = 1:nrow(data$x_test),
    depth_node = 0,
    node_var = NA,
    node_cutpoint_index = NA,
    left = NA,
    right = NA,
    parent_node = NA,
    ancestors = NA,
    terminal = TRUE,
    betas_vec = rep(0, length(unlist(data$basis_subindex)))
  )

  # Returning the node
  return(node)

}

# Get all the terminal nodes
get_terminals <- function(tree){

  # Return the name of the termianl nodes
  return(names(tree)[unlist(lapply(tree, function(x){x$terminal}),use.names =  TRUE)])
}

# Get nog terminal nodes
get_nogs <- function(tree){

  # Return the name of the termianl nodes
  non_terminal <- names(tree)[!unlist(lapply(tree, function(x){x$terminal}),use.names =  TRUE)]

  # In case there are non nonterminal nondes
  if(length(non_terminal)==0){
    return(non_terminal)
  }

  bool_nog <- vector("logical",length = length(non_terminal))
  for(i in 1:length(bool_nog)){
    # Checking if both children are terminal
    if( tree[[tree[[non_terminal[i]]]$left]]$terminal & tree[[tree[[non_terminal[i]]]$right]]$terminal) {
      bool_nog[i] <- TRUE
    }
  }

  return(  non_terminal[bool_nog])
}

# Getting the maximum node index number
get_max_node <- function(tree){

  # Return the name of the termianl nodes
  return(max(unlist(lapply(tree, function(x){x$node_number}),use.names =  TRUE)))
}




# A function to calculate the loglikelihood
nodeLogLike <- function(curr_part_res,
                        j_,
                        index_node,
                        data){

  # Subsetting the residuals
  curr_part_res_leaf <- curr_part_res[index_node]

  # Getting the number of observationsin the terminal node
  n_leaf <- length(index_node)
  d_basis <- length(j_)
  ones <- matrix(1,nrow = n_leaf)


  if(length(j_)==0){
    stop(" Node Log-likelihood: No variables")
  }

  # Using the Andrew's approach I would have
  mean_aux <- rep(0,length(curr_part_res_leaf))
  cov_aux <- matrix(0,nrow = length(index_node),ncol = length(index_node))
  # diag_tau_beta_inv <- diag(x = 1/unique(data$tau_beta), nrow = )

  for(jj in 1:length(j_)){
    # if(data$dif_order==0){
    #   stop("Never")
    #   # P_aux <- data$P[D_subset_index,D_subset_index]
    #   # for(jj in 1:length(j_)){
    #   #   P_aux[data$basis_subindex[[jj]],data$basis_subindex[[jj]]] <- data$tau_beta[jj]*data$P[data$basis_subindex[[jj]],data$basis_subindex[[jj]]]
    #   # }
    #   cov_aux <- diag(x = (data$tau^(-1)),nrow = n_leaf) + D_leaf%*%tcrossprod(diag_tau_beta_inv,D_leaf)
    # } else {

      # Adding the quantities with respect to the interaction
      if(j_[jj] <= length(data$dummy_x$continuousVars)){
        cov_aux <- cov_aux + data$B_train[[j_[jj]]][index_node,,drop = FALSE]%*%solve(data$P,t(data$B_train[[j_[jj]]][index_node,,drop = FALSE]))
      } else {
        cov_aux <- cov_aux + data$B_train[[j_[jj]]][index_node,,drop = FALSE]%*%solve(data$P_interaction,t(data$B_train[[j_[jj]]][index_node,,drop = FALSE]))
      }
    # }
  }
  cov_aux <- diag(x = (data$tau^(-1)),nrow = n_leaf) + cov_aux

  result <- mvnfast::dmvn(X = curr_part_res_leaf,mu = mean_aux,
                          sigma = cov_aux ,log = TRUE)


  return(c(result))

}



# Grow a tree
grow <- function(tree,
                 curr_part_res,
                 data){

  # Getting the maximum index number
  max_index <- get_max_node(tree)

  # Sampling a terminal node
  terminal_nodes <- get_terminals(tree)
  n_t_nodes <- length(terminal_nodes)
  nog_nodes <- get_nogs(tree)
  n_nog_nodes <- length(nog_nodes)
  g_node_name <- sample(terminal_nodes,size = 1)
  g_node <- tree[[g_node_name]]


  valid_terminal_node <- TRUE
  valid_count <- 0

  # acceptance_grid <- numeric(100)
  # for(kk in 1:100){
  while(valid_terminal_node){
    # Convinience while to avoid terminal nodes of 2

    # Sample a split var
    # ===== Uncomment this line below after ========
    p_var <- sample(1:NCOL(data$x_train),size = 1)
    # ==============================================
    # p_var <- 7

    # Selecting an available cutpoint from this terminal node
    valid_range_grow <- range(data$x_train[g_node$train_index,p_var])

    # Case of invalid range
    if(length(valid_range_grow)==0){
      return(tree)
    }

    # Subsetting the indexes of
    valid_cutpoint <- which(data$xcut_m[,p_var]>valid_range_grow[1] & data$xcut_m[,p_var]<valid_range_grow[2])

    # When there's no valid cutpoint on the sampled terminal node
    if(length(valid_cutpoint)==0){
      return(tree)
    }

    # Getting which cutpoints are valid and sample onde index
    sample_cutpoint <- sample(valid_cutpoint,
                              size = 1)
    # sample_cutpoint <- valid_cutpoint[kk]

    # Getting the left & right index
    left_index  <- data$all_var_splits[[p_var]][[sample_cutpoint]]$left_train[data$all_var_splits[[p_var]][[sample_cutpoint]]$left_train %in% g_node$train_index]
    right_index  <- data$all_var_splits[[p_var]][[sample_cutpoint]]$right_train[data$all_var_splits[[p_var]][[sample_cutpoint]]$right_train %in% g_node$train_index]

    left_test_index  <- data$all_var_splits[[p_var]][[sample_cutpoint]]$left_test[data$all_var_splits[[p_var]][[sample_cutpoint]]$left_test %in% g_node$test_index]
    right_test_index  <- data$all_var_splits[[p_var]][[sample_cutpoint]]$right_test[data$all_var_splits[[p_var]][[sample_cutpoint]]$right_test %in% g_node$test_index]



    # Verifying that the correct number was used
    if((length(left_index)+length(right_index))!=length(g_node$train_index)){
      stop("Something went wrong here --- train grown index doest match")
    }

    if((length(left_test_index)+length(right_test_index))!=length(g_node$test_index)){
      stop("Something went wrong here --- test grown index doest match")
    }


    # === Uncomment those lines after
    if( (length(left_index) > data$node_min_size) & (length(right_index)>data$node_min_size)){
      # Getting out of the while
      break
    } else {

      # Adding one to the counter
      valid_count = valid_count + 1

      # Stop trying to search for a valid cutpoint
      if(valid_count > 2) {
        valid_terminal_node = FALSE
        return(tree)
      }
    }
  }

  # For convinience we are going to avoid terminal nodes less than 2
  if( (length(left_index)<2) || (length(right_index) < 2)) {
    stop("Error of invalid terminal node")
  }

  # Calculating loglikelihood for the grown node, the left and the right node
  # Recover the g_node index

  if(!any(is.na(g_node$inter))){
    # node_index_var <- c(g_node$j,which( names(data$basis_subindex) %in% paste0(g_node$j,sort(g_node$inter))))
    node_index_var <- g_node$pred_vars # Now this information is stored here
  } else {
    # node_index_var <- g_node$j
    node_index_var <- g_node$pred_vars # Now this information is stored in that variable
  }

  g_loglike <- nodeLogLike(curr_part_res = curr_part_res,
                           j_ = node_index_var,
                           index_node = g_node$train_index,
                           data = data)


  left_loglike <-  nodeLogLike(curr_part_res = curr_part_res,
                               j_ = node_index_var,
                               index_node = left_index,
                               data = data)

  right_loglike <-  nodeLogLike(curr_part_res = curr_part_res,
                                j_ = node_index_var,
                                index_node = right_index,
                                data = data)

  # Calculating the prior
  prior_loglike <- log(data$alpha*(1+g_node$depth_node)^(-data$beta)) + # Prior of the grown node becoming nonterminal
    2*log(1-data$alpha*(1+g_node$depth_node+1)^(-data$beta)) - # plus the prior of the two following nodes being terminal
    log(1-data$alpha*(1+g_node$depth_node)^(-data$beta)) # minus the probability of the grown node being terminal

  # Transition prob
  log_trasition_prob  = log(0.3/(n_nog_nodes+1))-log(0.3/n_t_nodes)

  # Calculating the acceptance probability
  acceptance <- exp(-g_loglike+left_loglike+right_loglike+prior_loglike+log_trasition_prob)
  # acceptance_grid[kk] <- acceptance



  # par(mfrow=c(1,2))
  # plot(data$xcut_m[,2],(acceptance_grid), main = "Acceptance to split on X2", xlab = "X2", ylab = "Prob. Acceptance")

  if(data$stump) {
    acceptance <- acceptance*(-1)
  }

  # Getting the training the left and the right index for the the grown node
  if(stats::runif(n = 1)<acceptance){

    if(any(is.na(g_node$ancestors))){
      new_ancestors <- p_var
    } else {
      new_ancestors <- c(g_node$ancestors,p_var)
    }

    left_node <- list(node_number = max_index+1,
                      j = g_node$j,
                      pred_vars = g_node$pred_vars,
                      inter = g_node$inter,
                      isRoot = FALSE,
                      train_index = left_index,
                      test_index = left_test_index,
                      depth_node = g_node$depth_node+1,
                      node_var = p_var,
                      node_cutpoint_index = sample_cutpoint,
                      left = NA,
                      right = NA,
                      parent_node = g_node_name,
                      ancestors = new_ancestors,
                      terminal = TRUE,
                      betas_vec = g_node$betas_vec)

    right_node <- list(node_number = max_index+2,
                       j = g_node$j,
                       pred_vars = g_node$pred_vars,
                       inter = g_node$inter,
                       isRoot = FALSE,
                       train_index = right_index,
                       test_index = right_test_index,
                       depth_node = g_node$depth_node+1,
                       node_var = p_var,
                       node_cutpoint_index = sample_cutpoint,
                       left = NA,
                       right = NA,
                       parent_node = g_node_name,
                       ancestors = new_ancestors,
                       terminal = TRUE,
                       betas_vec = g_node$betas_vec)

    # Modifying the current node
    tree[[g_node_name]]$left = paste0("node",max_index+1)
    tree[[g_node_name]]$right = paste0("node",max_index+2)
    tree[[g_node_name]]$terminal = FALSE

    tree[[paste0("node",max_index+1)]] <- left_node
    tree[[paste0("node",max_index+2)]] <- right_node


  } else {

    # Do nothing

  }

  # Return the new tree
  return(tree)
}


# Adding interaction
add_interaction <- function(tree,
                            curr_part_res,
                            data){

  # Getting the maximum index number
  max_index <- get_max_node(tree)

  # Sampling a terminal node
  terminal_nodes <- get_terminals(tree)
  n_t_nodes <- length(terminal_nodes)
  nog_nodes <- get_nogs(tree)
  n_nog_nodes <- length(nog_nodes)
  g_node_name <- sample(terminal_nodes,size = 1)
  # g_node_name <- "node2"
  g_node <- tree[[g_node_name]]


  valid_terminal_node <- TRUE
  valid_count <- 0

  # Maybe need to change this when we start to work with continuous variables
  interaction_candidates <- (1:NCOL(data$x_train))[-g_node$j]

  # In case of the main effects are not included
  if(!(g_node$j %in% g_node$pred_vars)){
    interaction_candidates <- (1:NCOL(data$x_train))
  }

  while(length(interaction_candidates)!=0){

    # Sample a split var
    # ===== Uncomment this line below after ========
    p_var <- sample(interaction_candidates,size = 1)
    # p_var <- 10
    # ==============================================

    # In case of sampling the main effect
    if(p_var==g_node$j){
      int_name_ <- g_node$j
    } else {
      # Getting the interaction name
      int_name_ <- paste0(sort(c(p_var,g_node$j)),collapse = '')
    }

    # Making sure that not selecting a interaction that's already in the terminal node
    if(any(is.na(g_node$inter))){
      break
    } else {
      if(p_var %in% g_node$inter){
        index_candidate <- which(interaction_candidates %in% p_var)
        interaction_candidates <- interaction_candidates[-index_candidate]
      } else {
        break # Considering the case that the variable isnt there
      }
    }

  }

  # All interactions were already included.
  if(length(interaction_candidates)==0){
    return(tree)
  }

  # Calculating loglikelihood for the grown node, the left and the right node

  g_loglike <- nodeLogLike(curr_part_res = curr_part_res,
                           j_= g_node$pred_vars, # Here j_ refers to the main effect from that tree
                           index_node = g_node$train_index,
                           data = data)

  # Adding the interaction term on the list
  new_j <- sort(unique(c(g_node$pred_vars,which(names(data$basis_subindex) %in% int_name_))))

  # new_j <- c(10)
  new_loglike <-  nodeLogLike(curr_part_res = curr_part_res,
                              j_ = new_j,
                              index_node = g_node$train_index,
                              data = data)



  # Calculating the acceptance probability
  acceptance <- exp(-g_loglike+new_loglike)
  # acceptance_grid[kk] <- acceptance



  # par(mfrow=c(1,2))
  # plot(data$xcut_m[,2],(acceptance_grid), main = "Acceptance to split on X2", xlab = "X2", ylab = "Prob. Acceptance")

  if(data$stump) {
    acceptance <- acceptance*(-1)
  }

  if(g_node$j %in% interaction_candidates){
    Stop("Cannot have the main effect variable into the candidates of the model")
  }

  # Getting the training the left and the right index for the the grown node
  if(stats::runif(n = 1)<acceptance){

    # Modifying the current node (case is the first interaction)
    if(any(is.na(tree[[g_node_name]]$inter))){
      tree[[g_node_name]]$inter <- p_var
      tree[[g_node_name]]$pred_vars <- new_j
    } else {
      tree[[g_node_name]]$inter = sort(unique(c(tree[[g_node_name]]$inter,p_var))) # Case of previous interactions
      tree[[g_node_name]]$inter <- tree[[g_node_name]]$inter[!(tree[[g_node_name]]$inter %in% tree[[g_node_name]]$j)]# Excluding the main effect
      tree[[g_node_name]]$pred_vars <- new_j
    }

  } else {
    # Do nothing

  }

  # Return the new tree
  return(tree)
}


# Adding interaction
add_variable<- function(tree,
                            curr_part_res,
                            data){

  # Getting the maximum index number
  max_index <- get_max_node(tree)

  # Sampling a terminal node
  terminal_nodes <- get_terminals(tree)
  n_t_nodes <- length(terminal_nodes)
  nog_nodes <- get_nogs(tree)
  n_nog_nodes <- length(nog_nodes)
  g_node_name <- sample(terminal_nodes,size = 1)
  # g_node_name <- "node2"
  g_node <- tree[[g_node_name]]


  valid_terminal_node <- TRUE
  valid_count <- 0

  # Maybe need to change this when we start to work with continuous variables
  candidates <- (1:(length(data$basis_subindex)))[!(1:(length(data$basis_subindex)) %in% g_node$pred_vars)]


  # Sample a split var
  p_var <- sample(candidates,size = 1)

  # All interactions were already included.
  if(length(candidates)==0){
    return(tree)
  }

  # Calculating loglikelihood for the grown node, the left and the right node

  g_loglike <- nodeLogLike(curr_part_res = curr_part_res,
                           j_= g_node$pred_vars, # Here j_ refers to the main effect from that tree
                           index_node = g_node$train_index,
                           data = data)

  # Adding the interaction term on the list
  new_j <- unique(sort(c(g_node$pred_vars,p_var)))

  # new_j <- c(10)
  new_loglike <-  nodeLogLike(curr_part_res = curr_part_res,
                              j_ = new_j,
                              index_node = g_node$train_index,
                              data = data)



  # Calculating the acceptance probability
  acceptance <- exp(-g_loglike+new_loglike)
  # acceptance_grid[kk] <- acceptance



  # par(mfrow=c(1,2))
  # plot(data$xcut_m[,2],(acceptance_grid), main = "Acceptance to split on X2", xlab = "X2", ylab = "Prob. Acceptance")

  if(data$stump) {
    acceptance <- acceptance*(-1)
  }

  # Getting the training the left and the right index for the the grown node
  if(stats::runif(n = 1)<acceptance){

    # Modifying the current node (case is the first interaction)
    tree[[g_node_name]]$pred_vars <- new_j

  } else {
    # Do nothing

  }

  # Return the new tree
  return(tree)
}

# Pruning a tree
prune <- function(tree,
                  curr_part_res,
                  data){


  # Getting the maximum index number
  max_index <- get_max_node(tree)

  # Sampling a terminal node
  terminal_nodes <- get_terminals(tree)
  n_t_nodes <- length(terminal_nodes)
  nog_nodes <- get_nogs(tree)
  n_nog_nodes <- length(nog_nodes)

  # Just in case to avoid errors
  if(n_nog_nodes==0){
    return(tree)
  }

  # Selecting a node to be pruned
  p_node_name <- sample(nog_nodes,size = 1)
  p_node <- tree[[p_node_name]]

  # Getting the indexes from the left and right children from the pruned node
  children_left_index <- tree[[p_node$left]]$train_index
  children_right_index <- tree[[p_node$right]]$train_index
  children_left_ancestors <- tree[[p_node$left]]$ancestors
  children_right_ancestors <- tree[[p_node$right]]$ancestors

  # Calculating loglikelihood for the grown node, the left and the right node

  if(!any(is.na(p_node$inter))){
    # node_index_var <- c(p_node$j,which( names(data$basis_subindex) %in% paste0(p_node$j,sort(p_node$inter))))
    node_index_var <- p_node$pred_vars
  } else {
    # node_index_var <- p_node$j
    node_index_var <- p_node$pred_vars
  }

  p_loglike <- nodeLogLike(curr_part_res = curr_part_res,
                           index_node = p_node$train_index,
                           j_ = node_index_var,
                           data = data)


  p_left_loglike <-  nodeLogLike(curr_part_res = curr_part_res,
                                 index_node =  children_left_index,
                                 j_ = node_index_var,
                                 data = data)

  p_right_loglike <-  nodeLogLike(curr_part_res = curr_part_res,
                                  index_node = children_right_index,
                                  j_ = node_index_var,
                                  data = data)

  # Calculating the prior
  prior_loglike <- log(1-data$alpha*(1+p_node$depth_node)^(-data$beta)) - # Prior of the new terminal node
    log(data$alpha*(1+p_node$depth_node)^(-data$beta)) - # Prior of the grown node becoming nonterminal
    2*log(1-data$alpha*(1+p_node$depth_node+1)^(-data$beta))  # plus the prior of the two following nodes being terminal
  # minus the probability of the grown node being terminal

  # Transition prob
  log_trasition_prob  = log(0.3/(n_t_nodes))-log(0.3/n_nog_nodes)

  # Calculating the acceptance probability
  acceptance <- exp(p_loglike-p_left_loglike-p_right_loglike+prior_loglike+log_trasition_prob)

  # Getting the training the left and the right index for the the grown node
  if(stats::runif(n = 1)<acceptance){

    # Erasing the terminal nodes
    tree[[p_node$left]] <- NULL
    tree[[p_node$right]] <- NULL

    # Modifying back the pruned node
    tree[[p_node_name]]$left <- NA
    tree[[p_node_name]]$right <- NA
    tree[[p_node_name]]$terminal <- TRUE

  } else {
    # Do nothing
  }

  # Return the new tree
  return(tree)

}

# Pruning a tree
prune_interaction <- function(tree,
                              curr_part_res,
                              data){



  # Getting the maximum index number
  max_index <- get_max_node(tree)

  # Sampling a terminal node
  terminal_nodes <- get_terminals(tree)
  n_t_nodes <- length(terminal_nodes)

  t_with_inter <- names(which(sapply(terminal_nodes,function(node){(!all(is.na(tree[[node]]$inter))) & (length(tree[[node]]$node_vars)>1)})))


  if(length(t_with_inter)==0){
    # to avoid to waste an interaction
    std_prune <- prune(tree = tree,
                       curr_part_res = curr_part_res,data = data)
    return(std_prune)
  }

  # Selecting a node to be pruned
  p_node_name <- sample(t_with_inter,size = 1)
  p_node <- tree[[p_node_name]]


  # Calculating loglikelihood for the grown node, the left and the right node
  if(!any(is.na(p_node$inter))){
    # node_index_var <- c(p_node$j,which( names(data$basis_subindex) %in% paste0(sort(c(p_node$j,p_node$inter)),collapse = "")))
    node_index_var <- p_node$pred_vars
    inter_index_ <- p_node$inter

    # Sampling the new interactions subset
    if(length(inter_index_)==1){
      p_inter_index <- sample(p_node$pred_vars,size = 1)
      # p_inter_index <- p_node$inter
      if(p_inter_index %in% 1:NCOL(data$x_train)){
        new_p_inter <- p_node$inter
        new_node_index_var <- p_node$pred_vars[!(p_node$pred_vars %in% p_node$j)]
      } else {
        new_p_inter <- NA
        new_node_index_var <- p_node$pred_vars[p_node$pred_vars %in% p_node$j] # Getting only main effects
      }

      if(length(new_node_index_var)==0){
        stop("Error on the prune interaction")
      }
      # new_node_index_var <- c(p_node$j,which( names(data$basis_subindex) %in% paste0(p_node$j,sort(new_p_inter))))

    } else  {
      new_p_inter_index <- sort(sample(c(p_node$j,p_node$inter),size = 1,replace = FALSE))

      if(!(new_p_inter_index %in% p_node$j)){
        new_p_inter <- p_node$inter[!(p_node$inter %in% new_p_inter_index)]
        if(p_node$j %in% p_node$pred_vars){ # Checking if the main effect gonna be present in the predictors
          new_node_index_var <- c(p_node$j,which( names(data$basis_subindex) %in% paste0(sort(c(p_node$j,new_p_inter)),collapse = "")))
        } else {
          new_node_index_var <- which( names(data$basis_subindex) %in% paste0(sort(c(p_node$j,new_p_inter)),collapse = ""))
        }
      } else {
        new_p_inter <- p_node$inter
        new_node_index_var <- sort(unique(which( names(data$basis_subindex) %in% apply(sapply(new_p_inter,function(x){sort(c(p_node$j,x))}),2,function(y){paste0(y,collapse = "")}) )) )
      }
    }
  } else {
    stop('Prune interaction was called where there is no interaction')
  }

  p_loglike <- nodeLogLike(curr_part_res = curr_part_res,
                           j_ = p_node$pred_vars,
                           index_node = p_node$train_index,
                           data = data)


  new_p_loglike <-  nodeLogLike(curr_part_res = curr_part_res,
                                j_ = new_node_index_var,
                                index_node = p_node$train_index,
                                data = data)


  # Calculating the acceptance probability
  acceptance <- exp(-p_loglike+new_p_loglike)

  # Getting the training the left and the right index for the the grown node
  if(stats::runif(n = 1)<acceptance){

    # Erasing the terminal nodes
    tree[[p_node_name]]$inter <- new_p_inter
    tree[[p_node_name]]$pred_vars <- new_node_index_var

  } else {
    # Do nothing
  }

  # Return the new tree
  return(tree)

}


# Pruning a tree
prune_var <- function(tree,
                      curr_part_res,
                              data){



  # Getting the maximum index number
  max_index <- get_max_node(tree)

  # Sampling a terminal node
  terminal_nodes <- get_terminals(tree)
  n_t_nodes <- length(terminal_nodes)

  t_with_inter <- names(which(sapply(terminal_nodes,function(node){(length(tree[[node]]$pred_vars)>1)})))


  if(length(t_with_inter)==0){
    # to avoid to waste an interaction
    std_prune <- prune(tree = tree,
                       curr_part_res = curr_part_res,data = data)
    return(std_prune)
  }

  # Selecting a node to be pruned
  p_node_name <- sample(t_with_inter,size = 1)
  p_node <- tree[[p_node_name]]


  # Calculating loglikelihood for the grown node, the left and the right node
  p_index_var <- sample(1:length(p_node$pred_vars),size = 1)
  new_node_index_var <- sort(unique(p_node$pred_vars[-p_index_var]))


  p_loglike <- nodeLogLike(curr_part_res = curr_part_res,
                           j_ = p_node$pred_vars,
                           index_node = p_node$train_index,
                           data = data)


  new_p_loglike <-  nodeLogLike(curr_part_res = curr_part_res,
                                j_ = new_node_index_var,
                                index_node = p_node$train_index,
                                data = data)


  # Calculating the acceptance probability
  acceptance <- exp(-p_loglike+new_p_loglike)

  # Getting the training the left and the right index for the the grown node
  if(stats::runif(n = 1)<acceptance){

    # Erasing the terminal nodes
    tree[[p_node_name]]$pred_vars <- new_node_index_var

  } else {
    # Do nothing
  }

  # Return the new tree
  return(tree)

}

change_stump <- function(tree = tree,
                         curr_part_res = curr_part_res,
                         data = data){

  # Getting the stump
  c_node <- tree$node0

  # Proposing a change to the stump
  change_candidates <- which(!(1:NCOL(data$x_train) %in% c_node$pred_vars))

  # In case there's other proposal trees (only for 1-d case)
  if(length(change_candidates)==0){
    # Gettina grown tree
    grown_tree <- grow(tree = tree,
                       curr_part_res = curr_part_res,
                       data = data)
    return(grown_tree)
  }

  # Print
  new_ancestor <- sample(change_candidates,size = 1)


  stump_loglikelihood <- nodeLogLike(curr_part_res = curr_part_res,
                                     j_ = c_node$j,
                                     index_node = c_node$train_index,
                                     data = data)

  new_stump_loglikelihood <- nodeLogLike(curr_part_res = curr_part_res,
                                         j_ = new_ancestor,
                                         index_node = c_node$train_index,
                                         data = data)

  acceptance <- exp(-stump_loglikelihood+new_stump_loglikelihood)

  # Update the stump in case of changing the main var
  if(stats::runif(n = 1)<acceptance){
    # Modifying the node0
    tree$node0$j <- new_ancestor
    tree$node0$pred_vars <- new_ancestor
  }


  # Returning the new tree
  return(tree)


}


# Change a tree
change <- function(tree,
                   curr_part_res,
                   data){

  # Changing the stump
  if(length(tree)==1){
    change_stump_obj <- change_stump(tree = tree,
                                     curr_part_res = curr_part_res,
                                     data = data)
    return(change_stump_obj)
  }

  # Sampling a terminal node
  nog_nodes <- get_nogs(tree)
  n_nog_nodes <- length(nog_nodes)
  c_node_name <- sample(nog_nodes,size = 1)
  c_node <- tree[[c_node_name]]


  valid_terminal_node <- TRUE
  valid_count <- 0


  while(valid_terminal_node){
    # Convinience while to avoid terminal nodes of 2
    # Sample a split var
    p_var <- sample(1:ncol(data$x_train),size = 1)

    # Selecting an available cutpoint from this terminal node
    valid_range_grow <- range(data$x_train[c_node$train_index,p_var])

    # Subsetting the indexes of
    valid_cutpoint <- which(data$xcut_m[,p_var]>valid_range_grow[1] & data$xcut_m[,p_var]<valid_range_grow[2])

    # When there's no valid cutpoint on the sampled terminal node
    if(length(valid_cutpoint)==0){
      return(tree)
    }

    # Getting which cutpoints are valid and sample onde index
    sample_cutpoint <- sample(valid_cutpoint,
                              size = 1)

    # Getting the left & right index
    left_index  <- data$all_var_splits[[p_var]][[sample_cutpoint]]$left_train[data$all_var_splits[[p_var]][[sample_cutpoint]]$left_train %in% c_node$train_index]
    right_index  <- data$all_var_splits[[p_var]][[sample_cutpoint]]$right_train[data$all_var_splits[[p_var]][[sample_cutpoint]]$right_train %in% c_node$train_index]

    left_test_index  <- data$all_var_splits[[p_var]][[sample_cutpoint]]$left_test[data$all_var_splits[[p_var]][[sample_cutpoint]]$left_test %in% c_node$test_index]
    right_test_index  <- data$all_var_splits[[p_var]][[sample_cutpoint]]$right_test[data$all_var_splits[[p_var]][[sample_cutpoint]]$right_test %in% c_node$test_index]



    # Verifying that the correct number was used
    if((length(left_index)+length(right_index))!=length(c_node$train_index)){
      stop("Something went wrong here --- train grown index doest match")
    }

    if((length(left_test_index)+length(right_test_index))!=length(c_node$test_index)){
      stop("Something went wrong here --- test grown index doest match")
    }

    # Avoiding having terminal nodes with just one observation
    if( (length(left_index) > data$node_min_size) & (length(right_index)>data$node_min_size)){
      # Getting out of the while
      break
    } else {

      # Adding one to the counter
      valid_count = valid_count + 1

      # Stop trying to search for a valid cutpoint
      if(valid_count > 2) {
        valid_terminal_node = FALSE
        return(tree)
      }
    }
  }

  # For convinience we are going to avoid terminal nodes less than 2
  if( (length(left_index)<2) || (length(right_index) < 2)) {
    stop("Error of invalid terminal node")
  }


  # Getting the node_index var
  if(!any(is.na(c_node$inter))){
    # node_index_var <- c(c_node$j,which( names(data$basis_subindex) %in% paste0(c_node$j,sort(c_node$inter))))
    node_index_var <- c_node$pred_vars
  } else {
    node_index_var <- c_node$pred_vars
  }

  # Calculating loglikelihood for the new changed nodes and the old ones
  c_loglike_left <- nodeLogLike(curr_part_res = curr_part_res,
                                index_node = tree[[c_node$left]]$train_index,
                                j_ = node_index_var,
                                data = data)


  c_loglike_right <-  nodeLogLike(curr_part_res = curr_part_res,
                                  index_node = tree[[c_node$right]]$train_index,
                                  j_ =  node_index_var,
                                  data = data)

  # Calculating a new ancestors left and right
  old_p_var <- tree[[c_node$left]]$node_var

  # Storing new left and right ancestors
  new_left_ancestors <- tree[[c_node$left]]$ancestors
  new_left_ancestors[length(new_left_ancestors)] <- p_var

  new_right_ancestors <- tree[[c_node$right]]$ancestors
  new_right_ancestors[length(new_right_ancestors)] <- p_var


  new_c_loglike_left <-  nodeLogLike(curr_part_res = curr_part_res,
                                     index_node = left_index,
                                     j = node_index_var,
                                     data = data)

  new_c_loglike_right <-  nodeLogLike(curr_part_res = curr_part_res,
                                      index_node = right_index,
                                      j =  node_index_var,
                                      data = data)


  # Calculating the acceptance probability
  acceptance <- exp(new_c_loglike_left+new_c_loglike_right-c_loglike_left-c_loglike_right)

  # Getting the training the left and the right index for the the grown node
  if(stats::runif(n = 1,min = 0,max = 1)<acceptance){

    # Updating the left and the right node
    # === Left =====
    tree[[c_node$left]]$node_var <- p_var
    tree[[c_node$left]]$node_cutpoint_index <- sample_cutpoint
    tree[[c_node$left]]$train_index <- left_index
    tree[[c_node$left]]$test_index <- left_test_index
    tree[[c_node$left]]$ancestors <- new_left_ancestors

    #==== Right ====
    tree[[c_node$right]]$node_var <- p_var
    tree[[c_node$right]]$node_cutpoint_index <- sample_cutpoint
    tree[[c_node$right]]$train_index <- right_index
    tree[[c_node$right]]$test_index <- right_test_index
    tree[[c_node$right]]$ancestors <- new_right_ancestors

  } else {
    # Do nothing
  }

  # Return the new tree
  return(tree)

}

# Change interaction
change_interaction <-  function(tree,
                                curr_part_res,
                                data){


  # Getting the maximum index number
  max_index <- get_max_node(tree)

  # Sampling a terminal node
  terminal_nodes <- get_terminals(tree)
  n_t_nodes <- length(terminal_nodes)

  t_with_inter <- names(which(sapply(terminal_nodes,function(node){!all(is.na(tree[[node]]$inter))})))

  # For the cases where there's no interactions
  if(length(t_with_inter)==0){
    std_change <- change(tree = tree,
                         curr_part_res = curr_part_res,
                         data = data)
    return(std_change)
  }

  # Selecting a node to be pruned
  c_node_name <- sample(t_with_inter,size = 1)
  c_node <- tree[[c_node_name]]


  # Calculating loglikelihood for the grown node, the left and the right node
  if(!any(is.na(c_node$inter))){
    node_index_var <- c(c_node$j,which( names(data$basis_subindex) %in% paste0(c_node$j,sort(c_node$inter))))
    inter_index_ <- c_node$inter

    # Sampling the new interactions subset
    if(length(inter_index_)==1){
      c_inter_index <- c_node$inter
      new_interaction_candidates <- (1:NCOL(data$x_train))[-c(c_node$j,c_node$inter)]

      if(length(new_interaction_candidates)==0){
        stop('Change interaction error')
      }

      new_c_inter <-   sample(new_interaction_candidates,size = 1)
      aux_names <- apply(sapply(new_c_inter,function(x){sort(c(c_node$j,x))}),2,function(y){paste0(y,collapse = "")})

      if(c_node$j %in% c_node$pred_vars)  {
        new_node_index_var <- c(c_node$j,which( names(data$basis_subindex) %in% aux_names)) #
      } else {
        new_node_index_var <- which( names(data$basis_subindex) %in% aux_names) #
      }

    } else  {
      inter_to_change <- sample(c_node$inter,size = 1)
      new_interaction_candidates <- which(!(1:NCOL(data$x_train)) %in% c(c_node$inter,c_node$j))
      if(length(new_interaction_candidates)==0){
        return(tree)
      }
      new_c_inter_single  <- sample(new_interaction_candidates,size = 1)
      inter_to_change_index <- which(c_node$inter %in% inter_to_change)
      new_c_inter <- c_node$inter
      new_c_inter[inter_to_change_index] <- new_c_inter_single
      new_c_inter <- sort(new_c_inter)
      aux_names <- apply(sapply(new_c_inter,function(x){sort(c(c_node$j,x))}),2,function(y){paste0(y,collapse = "")})
      if(c_node$j %in% c_node$pred_vars)  {
        new_node_index_var <- c(c_node$j,which( names(data$basis_subindex) %in% aux_names)) #
      } else {
        new_node_index_var <- which( names(data$basis_subindex) %in% aux_names) #
      }
    }

  } else {
    stop('Change was called where there is no interaction')
  }

  c_loglike <- nodeLogLike(curr_part_res = curr_part_res,
                           j_ = c_node$pred_vars,
                           index_node = c_node$train_index,
                           data = data)


  new_c_loglike <-  nodeLogLike(curr_part_res = curr_part_res,
                                j_ = new_node_index_var,
                                index_node = c_node$train_index,
                                data = data)


  # Calculating the acceptance probability
  acceptance <- exp(-c_loglike+new_c_loglike)

  # Getting the training the left and the right index for the the grown node
  if(stats::runif(n = 1)<acceptance){

    # Erasing the terminal nodes
    tree[[c_node_name]]$inter <- new_c_inter
    tree[[c_node_name]]$pred_vars <- new_node_index_var
  } else {
    # Do nothing
  }

  # Return the new tree
  return(tree)

}


# Change interaction
change_var <-  function(tree,
                                curr_part_res,
                                data){


  # Getting the maximum index number
  max_index <- get_max_node(tree)

  # Sampling a terminal node
  terminal_nodes <- get_terminals(tree)
  n_t_nodes <- length(terminal_nodes)

  t_with_inter <- names(which(sapply(terminal_nodes,function(node){length(tree[[node]]$pred_vars)>1})))

  # For the cases where there's no interactions
  if(length(t_with_inter)==0){
    std_change <- change(tree = tree,
                         curr_part_res = curr_part_res,
                         data = data)
    return(std_change)
  }

  # Selecting a node to be pruned
  c_node_name <- sample(t_with_inter,size = 1)
  c_node <- tree[[c_node_name]]


  # Choosing a new candidate
  c_node_new <- sample(which(!(1:length(data$basis_subindex)) %in% c_node$pred_vars),size = 1)

  new_node_index_var <- c_node$pred_vars
  new_node_index_var[sample(1:length(c_node$pred_vars),size = 1)] <- c_node_new
  new_node_index_var <- sort(new_node_index_var)

  c_loglike <- nodeLogLike(curr_part_res = curr_part_res,
                           j_ = c_node$pred_vars,
                           index_node = c_node$train_index,
                           data = data)


  new_c_loglike <-  nodeLogLike(curr_part_res = curr_part_res,
                                j_ = new_node_index_var,
                                index_node = c_node$train_index,
                                data = data)


  # Calculating the acceptance probability
  acceptance <- exp(-c_loglike+new_c_loglike)

  # Getting the training the left and the right index for the the grown node
  if(stats::runif(n = 1)<acceptance){
    # Erasing the terminal nodes
    tree[[c_node_name]]$pred_vars <- new_node_index_var
  } else {
    # Do nothing
  }

  # Return the new tree
  return(tree)

}

# ============
# Update Betas
# ============
updateBetas <- function(tree,
                        curr_part_res,
                        data,
                        trees_fit,
                        j){


  # Getting the terminals
  t_nodes_names <- get_terminals(tree)


  # Getting the current prediction for that tree
  y_hat_train <- matrix(0,nrow = nrow(data$x_train),ncol = length(data$basis_subindex))
  y_hat_test <- matrix(0,nrow = nrow(data$x_test),ncol = length(data$basis_subindex))

  for(i in 1:length(t_nodes_names)){


    # Select the current terminal node
    cu_t <- tree[[t_nodes_names[i]]]
    # THIS LINE IS COMPLETELY IMPORTANT BECAUSE DEFINE THE ANCESTEORS BY j only

    # Node VAR INDEX NOW IS INCLUDED IN THE TERMINAL NODE
    # ==============================
    # if(!any(is.na(cu_t$inter))){
    #   node_index_var <- c(cu_t$j,which( names(data$basis_subindex) %in% paste0(cu_t$j,sort(cu_t$inter))))
    # } else {
    #   node_index_var <- cu_t$j
    # }
    # ===========================
    # The lines above are summarised here
    node_index_var <- cu_t$pred_vars

    # Selecting the actually parameters subsetting
    basis_dim <- NCOL(data$P)
    basis_dim_interaction <- NCOL(data$P_interaction)
    n_leaf <- length(cu_t$train_index)
    diag_leaf <- diag(nrow = n_leaf)
    diag_basis <- diag(nrow = basis_dim)


    #  Calculating the quantities need to the posterior of \beta
    # == Starting to iterate over those coefficients ==========#
    for(jj in 1:length(node_index_var)){


      leaf_basis_subindex <- unlist(data$basis_subindex[node_index_var[jj]]) # Recall to the unique() here too

      # RES_LEAF also need to updated here from the new_curr_part_res
      old_betas <- matrix(tree[[t_nodes_names[i]]]$betas_vec[leaf_basis_subindex],nrow = 1)

      res_leaf <- matrix(curr_part_res[cu_t$train_index], ncol=1) - (trees_fit[j,cu_t$train_index] - tcrossprod(data$B_train[[node_index_var[jj]]][cu_t$train_index,,drop=FALSE],old_betas))

      # Getting the index for the vector of betas

      b_ <- crossprod(data$B_train[[node_index_var[jj]]][cu_t$train_index,,drop=FALSE],res_leaf)

      data_tau_beta_diag <- rep(data$tau_beta[node_index_var], NCOL(data$B_train[[node_index_var[jj]]])) # Don't really use this
      if(node_index_var[jj]<=length(data$dummy_x$continuousVars)){
        U_ <- data$P*data$tau_beta[node_index_var[jj]]
      } else {
        U_ <- data$P_interaction*data$tau_beta[node_index_var[jj]]
      }

      Q_ <- (crossprod(data$B_train[[node_index_var[jj]]]) + data$tau^(-1)*U_)
      Q_inv_ <- chol2inv(chol(Q_))
      # Q_inv_ <- solve(Q_)

      # Storing the old betas
      # See that I also creating a vector with the new betas
      new_betas <- mvnfast::rmvn(n = 1,mu = Q_inv_%*%b_,sigma = (data$tau^(-1))*Q_inv_)
      tree[[t_nodes_names[i]]]$betas_vec[leaf_basis_subindex] <- new_betas
      new_betas <- matrix(new_betas,nrow = 1)
      # Updating the residuals
      new_partial_pred <- tcrossprod(data$B_train[[node_index_var[jj]]][cu_t$train_index,,drop=FALSE],new_betas)
      # Need to update the trees fit!
      trees_fit[j,cu_t$train_index] <- trees_fit[j,cu_t$train_index] - tcrossprod(data$B_train[[node_index_var[jj]]][cu_t$train_index,,drop=FALSE],old_betas) + new_partial_pred

      y_hat_train[cu_t$train_index,node_index_var[jj]] <- new_partial_pred
      y_hat_test[cu_t$test_index,node_index_var[jj]] <- tcrossprod(data$B_test[[node_index_var[jj]]][cu_t$test_index,,drop=FALSE],new_betas)
    }

  }

  # Returning the tree
  return(list(tree = tree,
              y_hat_train = y_hat_train,
              y_hat_test = y_hat_test))

}


# =================
# Update \tau_betas
# =================
update_tau_betas_j <- function(forest,
                               data){


  # if(data$dif_order!=0){
  #   stop("Do not update tau_beta for peanalised version yet")
  # }

  # Setting some default hyperparameters
  # a_tau_beta <- d_tau_beta <- 0.1
  # Setting some default hyperparameters
  # a_tau_beta <- 0.1
  # d_tau_beta <- 0.1

  # Setting some default hyperparameters
  a_tau_beta <- data$a_tau_beta_j
  d_tau_beta <- data$d_tau_beta_j

  tau_b_shape <- 0.0
  tau_b_rate <- 0.0


  if(data$interaction_term){
    tau_b_shape <- numeric(NCOL(data$x_train)+NCOL(data$interaction_list))
    tau_b_rate <- numeric(NCOL(data$x_train)+NCOL(data$interaction_list))
    tau_beta_vec_aux <- numeric(NCOL(data$x_train)+NCOL(data$interaction_list))
  } else{
    tau_b_shape <- numeric(NCOL(data$x_train))
    tau_b_rate <- numeric(NCOL(data$x_train))
    tau_beta_vec_aux <- numeric(NCOL(data$x_train))
  }

  # Iterating over all trees
  for(i in 1:length(forest)){

    # Getting terminal nodes
    t_nodes_names <- get_terminals(forest[[i]])
    n_t_nodes <- length(t_nodes_names)

    # Iterating over the terminal nodes
    for(j in 1:length(t_nodes_names)){

      cu_t <- forest[[i]][[t_nodes_names[j]]]
      cu_t$ancestors <- cu_t$j
      # var_ <- cu_t$ancestors
      #
      #
      # # Getting the interactions as well
      # if(!any(is.na(cu_t$inter))){
      #   interaction_index <- cu_t$inter
      #   interaction_index <- sapply(interaction_index,function(x){sort(c(cu_t$j,x))})
      #   for(ii in 1:NCOL(interaction_index)){
      #     var_ <- c(var_,paste0(interaction_index[,ii],collapse = ""))
      #   }
      #   # for(var_ in 1:cu_t$ancestors){
      #   var_ <- which(names(data$basis_subindex) %in% var_)
      # }

      # All the information from var_ now is summarised inside the element from ht enode pred_vars
      var_ <- cu_t$pred_vars


      # Getting ht leaf basis
      for(kk in 1:length(var_)){
        leaf_basis_subindex <- unlist(data$basis_subindex[var_[kk]]) # Recall to the unique() function here
        p_ <- length(leaf_basis_subindex)
        betas_mat_ <- matrix(cu_t$betas_vec[leaf_basis_subindex],nrow = p_)
        tau_b_shape[var_[kk]] <- tau_b_shape[var_[kk]] + p_

        if(var_[kk] <= NCOL(data$x_train)){
          tau_b_rate[var_[kk]] <- tau_b_rate[var_[kk]] + c(crossprod(betas_mat_,crossprod(data$P,betas_mat_)))
        } else {
          tau_b_rate[var_[kk]] <- tau_b_rate[var_[kk]] + c(crossprod(betas_mat_,crossprod(data$P_interaction,betas_mat_)))

        }
      }
      # }

    }


  }

  if(data$interaction_term){
    for(j in 1:(NCOL(data$x_train)+NCOL(data$interaction_list)) ){
      tau_beta_vec_aux[j] <- rgamma(n = 1,
                                    shape = 0.5*tau_b_shape[j] + a_tau_beta,
                                    rate = 0.5*tau_b_rate[j] + d_tau_beta)

    }
  } else {
    for(j in 1:NCOL(data$x_train)){
      tau_beta_vec_aux[j] <- rgamma(n = 1,
                                    shape = 0.5*tau_b_shape[j] + a_tau_beta,
                                    rate = 0.5*tau_b_rate[j] + d_tau_beta)

    }
  }

  return(tau_beta_vec_aux)

}


update_tau_betas <- function(forest,
                             data){

  if(data$dif_order!=0){
    stop("Do not update tau_beta for peanalised version yet")
  }

  # Setting some default hyperparameters
  a_tau_beta <- d_tau_beta <- 0.1
  tau_b_shape <- 0.0
  tau_b_rate <- 0.0


  # Iterating over all trees
  for(i in 1:length(forest)){

    # Getting terminal nodes
    t_nodes_names <- get_terminals(forest[[i]])
    n_t_nodes <- length(t_nodes_names)

    # Iterating over the terminal nodes
    for(j in 1:length(t_nodes_names)){

      cu_t <- forest[[i]][[t_nodes_names[j]]]
      leaf_basis_subindex <- unlist(data$basis_subindex[unique(cu_t$j)]) # Recall to the unique() function here

      if(!is.null(cu_t$betas_vec)){
        tau_b_shape <- tau_b_shape + length(leaf_basis_subindex)
        tau_b_rate <- tau_b_rate + c(crossprod(cu_t$betas_vec[leaf_basis_subindex]))
      }

    }


    tau_beta_vec_aux <- rgamma(n = 1,
                               shape = 0.5*tau_b_shape + a_tau_beta,
                               rate = 0.5*tau_b_rate + d_tau_beta)
  }

  return(tau_beta_vec_aux)

}


# ===================
# Updating the \delta
# ===================

# A function to get predictions
getPredictions <- function(tree,
                           data){

  # Creating the vector to hold the values of the prediction
  if(data$interaction_term){
    y_hat <- matrix(0, nrow = nrow(data$x_train), ncol = NCOL(data$x_train)+NCOL(data$interaction_list))
    y_hat_test <- matrix(0,nrow(data$x_test), ncol = NCOL(data$x_test)+NCOL(data$interaction_list))
  } else {
    y_hat <- matrix(0, nrow = nrow(data$x_train), ncol = ncol(data$x_train))
    y_hat_test <- matrix(0,nrow(data$x_test), ncol = ncol(data$x_test))
  }

  # Getting terminal nodes
  t_nodes <- get_terminals(tree = tree)
  n_t_nodes <- length(t_nodes)

  for(i in 1:n_t_nodes){


    # Getting the current terminal node
    cu_t <- tree[[t_nodes[[i]]]]
    leaf_train_index <- cu_t$train_index
    leaf_test_index <- cu_t$test_index
    # leaf_ancestors <- unique(tree[[t_nodes[[i]]]]$ancestors) # recall the unique() argument here


    #===========
    # Don't need to get this information in this way, all is stored in the .$node_vars
    # if(!any(is.na(cu_t$inter))){
    #   node_index_var <- c(cu_t$j,which( names(data$basis_subindex) %in% paste0(cu_t$j,sort(cu_t$inter))))
    # } else {
    #   node_index_var <- cu_t$j
    # }
    #===========

    # Getting the variables used in the model
    node_index_var <- cu_t$pred_vars

    leaf_ancestors <- node_index_var # here isnt really the ancestors, but the variables that are being used

    leaf_basis_subindex <- data$basis_subindex[leaf_ancestors]

    # This test doesn't make sense anymore
    # # Test unit
    # if(length(leaf_ancestors)!=length(leaf_basis_subindex)){
    #   stop("Error on the getPredictions function")
    # }

    # Only add the marginal effects if the variables are within that terminal node
    if(length(leaf_basis_subindex)!=0){
      for(k in 1:length(leaf_basis_subindex)){

        y_hat[leaf_train_index,leaf_ancestors[k]] <- y_hat[leaf_train_index,leaf_ancestors[k]] + data$D_train[leaf_train_index,leaf_basis_subindex[[k]], drop = FALSE]%*%tree[[t_nodes[i]]]$betas_vec[leaf_basis_subindex[[k]]]
        y_hat_test[leaf_test_index,leaf_ancestors[k]] <- y_hat_test[leaf_test_index,leaf_ancestors[k]] + data$D_test[leaf_test_index,leaf_basis_subindex[[k]], drop = FALSE]%*%tree[[t_nodes[i]]]$betas_vec[leaf_basis_subindex[[k]]]

      }
    }

  }

  # Returning both training and test set predictions
  return(list(y_train_hat = y_hat,
              y_hat_test = y_hat_test))

}

# Updating tau
update_tau <- function(y_train_hat,
                       data){

  # Sampling a tau value
  n_ <- nrow(data$x_train)
  tau_sample <- stats::rgamma(n = 1,shape = 0.5*n_+data$a_tau,rate = 0.5*crossprod((data$y_train-y_train_hat))+data$d_tau)

  return(tau_sample)

}


