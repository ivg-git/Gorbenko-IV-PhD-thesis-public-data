select_top_n_nodes <- function(g,n){
  deg <-  degree(g)
  bet <- betweenness(g)
  clo <- closeness(g)
  le <- local_efficiency(g)
  hc <- harmonic_centrality(g)
  df <- tibble(deg = deg,
               tair = names(deg))
  db <- tibble(bet = bet,
               tair = names(bet))
  dc <- tibble(clo = clo,
               tair = names(clo))
  dl <- tibble(le = le,
               tair = names(le))
  dh <- tibble(hc = hc,
               tair = names(hc))
  df <- df[order(df$deg, decreasing = T),]
  db <- db[order(db$bet, decreasing = T),]
  dc <- dc[order(dc$clo, decreasing = T),]
  dl <- dl[order(dl$le, decreasing = T),]
  dh <- dh[order(dh$hc, decreasing = T),]
  df <- df[c(1:n),]
  db <- db[c(1:n),]
  dc <- dc[c(1:n),]
  dl <- dl[c(1:n),]
  dh <- dh[c(1:n),]
  df <- pivot_longer(df,cols = 1, names_to = "centrality", values_to = "value")
  db <- pivot_longer(db,cols = 1, names_to = "centrality", values_to = "value")
  dc <- pivot_longer(dc,cols = 1, names_to = "centrality", values_to = "value")
  dl <- pivot_longer(dl,cols = 1, names_to = "centrality", values_to = "value")
  dh <- pivot_longer(dh,cols = 1, names_to = "centrality", values_to = "value")
  result <- rbind(df,db,dc,dl,dh) %>% filter(value > 0) %>% pivot_wider(names_from = "centrality", values_from = "value")
  return(result)
}
