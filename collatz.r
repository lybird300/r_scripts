
collatz <- function(start_number){
	s <- start_number
	c <- 0
	all_s <- s
	while (s !=1) {
	ifelse(s%%2 == 0, s <- s/2, s <- 3*s+1)
	c <- c+1
	all_s <- append(all_s,s)
	}
	return(c,all_s)
}