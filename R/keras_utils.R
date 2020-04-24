KL_Annealing <- R6::R6Class(
  "KL_Annealing",
  inherit = KerasCallback,
  public = list(
    kl_weight = 0,
    anneal_iter = 1000,

    initialize = function(kl_weight, anneal_iter) {
      self$kl_weight <- kl_weight
      self$anneal_iter <- anneal_iter
    },
    on_batch_end = function(batch, logs = list()) {
      step <- 1/self$anneal_iter
      k_set_value(self$kl_weight, min(k_get_value(self$kl_weight) + step, 1))
    }
  )
)

create_activations <- function(nonlinear, act) {
  if (length(act) == 1) {
    act <- lapply(1:length(nonlinear), function(x) act)
  }
  act.list <- vector("list", length(nonlinear))
  for (i in seq_along(nonlinear)) {
    if (nonlinear[i]) {
      act.list[[i]] <- act[[i]]
    } else {
      act.list[[i]] <- layer_activation(activation = NULL)
    }
  }
  act.list
}

make_decoder <- function(dec, act,
                         bias, n_init, out_dim, name = "d", last_act = NULL,
                         include_last = TRUE) {
  decoder_layers <- list()
  for (i in seq_along(dec)) {
    decoder_layers[[i]] <- layer_dense(units = dec[i], activation = act[[i]],
                                       use_bias = bias, kernel_initializer = n_init,
                                       name = paste0(name, i))
  }
  if (include_last) {
    decoder_layers[[length(dec)+1]] <- layer_dense(units = out_dim, activation = last_act,
                                                   kernel_initializer = n_init,
                                                   use_bias = bias, name = paste0(name, "_out"))
  }
  return(decoder_layers)
}


make_bn <- function(dec, batch_norm2 = FALSE) {
  bn_layers <- list()
  for (i in seq_along(dec)) {
    bn_layers[[i]] <- layer_batch_normalization(center = batch_norm2, scale = batch_norm2,
                                                name = paste0("bd", i))
  }
  return(bn_layers)
}
