

init_little_groups <-function () {

  little_groups_2Oh <- vector ( mode = "list" )

  little_groups_2Oh[["2Oh"]]  <- list( 
                                    d=c( 0,  0,  0), 
                                    group="2Oh", 
                                    rid=c( 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48), 
                                    rmid=c( 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48), 
                                    irreps=c("A1g","A2g","Eg","T1g","T2g","G1g","G2g","Hg","A1u","A2u","Eu","T1u","T2u","G1u","G2u","Hu"), 
                                    irreps_dim=c(1,1,2,3,3,2,2,4,1,1,2,3,3,2,2,4) )

  little_groups_2Oh[["2C4v"]] <- list( 
                                    d=c( 0,  0,  1), 
                                    group="2C4v", 
                                    rid=c( 1,4,7,10,13,16,19,48), 
                                    rmid=c( 2,3,5,6,38,39,44,45), 
                                    irreps=c("A1","A2","B1","B2","E","G1","G2" ), 
                                    irreps_dim=c(1,1,1,1,2,2,2) )

  little_groups_2Oh[["2C2v"]] <- list( 
                                    d=c( 1,  1,  0), 
                                    group="2C2v", 
                                    rid=c( 1,38,44,48), 
                                    rmid=c( 4,7,39,45), 
                                    irreps=c( "A1","A2","B1","B2","G1" ), 
                                    irreps_dim=c(1,1,1,1,2) )

  little_groups_2Oh[["2C3v"]] <- list( 
                                    d=c( 1,  1,  1), 
                                    group="2C3v", 
                                    rid=c( 1,20,24,28,32,48), 
                                    rmid=c( 37,39,41,43,45,47), 
                                    irreps=c("A1","A2","K1","K2","E","G1" ), 
                                    irreps_dim=c(1,1,1,1,2,2) )

  return( little_groups_2Oh )
}  # end of init_little_groups
