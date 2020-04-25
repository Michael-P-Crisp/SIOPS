
   
module si_stats

implicit none

contains


!--------Get a variety of performance metrics for each investigation, typically some form of average relating to differential settlement --------- 
!invoked by process_si and process_si_multi
subroutine calc_perf(ninv,nrep_MC,costmetric,pilecost,testcost,inv_test,inv_depths,inv_bh,invcount,goodloc,failurevals,avesides,costs,diffset,fcost,pcost,icost,probfail,avediff,diffgeo,diffgeo2,si_performance)

		!input
		integer, intent(in) :: ninv,nrep_MC,costmetric !number of investigations, MC realisations, choice of objective function, 
		real, intent(in) :: pilecost,testcost(:)  !pile and borehole cost per metre
		integer, intent(in) :: inv_test(:,:),inv_depths(:,:,:),inv_bh(:) !test type, borehole depth and number of boreholes for each investigation
		integer, intent(in) :: goodloc(:,:), invcount(:) !locations and number of valid MC realisations for each investigation
		real, intent(in) :: costs(:,:),avesides(:,:),diffset(:,:) !failure costs, average pile length and differential settlements for each investigation and MC realisation
		real, intent(in) :: failurevals(:) !upper and lower thresholds of differential settlement, where the lower corresponds to failure.
		
		!output
		real, intent(out) :: fcost(:),pcost(:),icost(:) !average failure cost, pile construction cost, investigation cost for each investigation
		real, intent(out) :: probfail(:), avediff(:) !probability of failure, average differential settlement for each investigation
		real, intent(out) :: diffgeo(:),diffgeo2(:) !either a geometric mean or geometric standard deviation above the mean depending on costmetric
		real, intent(out) :: si_performance(:) !values of chosen objective function
        
        !real probfail2(size(probfail)) !calculate probability of failure using the lognormal statistics
		
		!local variables
		real tempvec(nrep_MC) !temporary array
		real tempmean !temp mean
		real tempsd !temp standard deviation
		logical nonzero(nrep_MC)
		integer num_nonzero
		real, parameter :: tol = epsilon(tempmean) !tiny(tempmean) * 10000 !tolerance threshold for near-zero differential settlements
		
		integer i,bh !loop counter


            !add average costs together for each investigation

            do i=1,ninv
                
                
    
                !calculate average failure cost
                fcost(i) = sum(costs(goodloc(:invcount(i),i),i),1)/invcount(i)

                !add in average pile cost
                pcost(i) = pilecost*sum(avesides(goodloc(:invcount(i),i),i),1)/invcount(i) !sides(1,:,:)
    
                !add in site investigation costs
                icost(i) = 0
                do bh = 1,inv_bh(i)
                    icost(i) = icost(i) + inv_depths(bh,i,2)*testcost(inv_test(bh,i))
                end do


                
                !Treat differential settlements as normally distributed; can get arithmetic average plus or minus standard deviations
                !tempvec(:invcount(i)) = diffset(goodloc(:invcount(i),i),i)
                !tempmean = sum(tempvec(:invcount(i)),1)/invcount(i)
                !tempvec(:invcount(i)) = tempvec(:invcount(i)) - tempmean
                !tempsd = sum(tempvec(:invcount(i))**2)
                !tempsd = sqrt(tempsd/invcount(i))
                tempmean = sum(diffset(goodloc(:invcount(i),i),i))/invcount(i)
                avediff(i) = tempmean !+ tempsd * costmetric
                !avediff2(i) = sum(diffset(goodloc(:invcount(i),i),i)**2)/invcount(i)     
                
            	!Treat differential settlements as LOGnormally distributed; can get geometric average plus or minus geometric standard deviations
            	!zero and near-zero values are removed because the final result becomes 0 (actually, undefined here) and does not reflect the variabiilty

                tempvec(:invcount(i)) = diffset(goodloc(:invcount(i),i),i)			!extract valid values
                nonzero(:invcount(i)) = tempvec(:invcount(i)) > tol	!find where values area greater than 0
                num_nonzero = count(nonzero(:invcount(i)))							!count non-zero values
                if(num_nonzero == 0) then											!if there are only zero values, set to zero
                	diffgeo(i) = 0
                    diffgeo2(i) = 0
                else
                    tempvec(:num_nonzero) = pack(tempvec(:invcount(i)),nonzero(:invcount(i)))	!extract non-zero values
                    tempvec(:num_nonzero) = log(tempvec(:num_nonzero))							!take log for geometric calculations
                	tempmean = sum(tempvec(:num_nonzero),1)/num_nonzero
                	tempvec(:invcount(i)) = tempvec(:num_nonzero) - tempmean
                	tempsd = sum(tempvec(:num_nonzero)**2)
                	tempsd = sqrt(tempsd/num_nonzero)
                	!diffgeo0(i) = tempmean
                	diffgeo(i) = exp(tempmean) * exp(tempsd) ** costmetric
                    diffgeo2(i) = exp(tempmean)
                	!diffgeo2(i) = tempmean * tempsd ** 2
				end if

!                 if(exp(tempmean) <= epsilon(tempmean)) then !if the geometric mean is zero (which occurs if any value is zero), set it explicitly, or else the SD is a NaN
! 					diffgeo(i) = 0
! 				else
! 					tempvec(:invcount(i)) = tempvec(:invcount(i)) - tempmean
!                 	tempsd = sum(tempvec(:invcount(i))**2)
!                 	tempsd = sqrt(tempsd/invcount(i))
!                 	!diffgeo0(i) = tempmean
!                 	diffgeo(i) = exp(tempmean) * exp(tempsd) ** costmetric
!                 	!diffgeo2(i) = tempmean * tempsd ** 2
! 				end if
       
                !probability of failure, the proportion of differential settlements above the failure threshold
                probfail(i) = 100*real(count(diffset(goodloc(:invcount(i),i),i) > failurevals(1)))/invcount(i)
            
                
       			!---A potentially more precise calculation of probability of failure, particularly with a low number of MC realisations, or at the tail end of the distribution
       			!---assumes that the differential settlement is lognormally-distributed
                !---- testing has shown it to be too inaccurate
       			!probfail2(i) = 1 - logcdf(failurevals(1),tempmean,tempsd)
       
       
                !evalsmean(i) = sum(evals(goodloc(:invcount(i),i),i))/invcount(i)
                !evalsmean(i) = sum(evals(:,i))/nrep_MC

            end do
            

            
                !Choose appropriate performance metric for the EA
				if (costmetric == -4) then
					si_performance = fcost      !use cost information for site investigation performance
	
				else if (costmetric == -3) then
					si_performance = fcost + pcost + icost      !use cost information for site investigation performance
                    
                else if (costmetric == -2) then
					si_performance = avediff    !average weighted differential settlemen
		
				else if (costmetric == -1) then 
					si_performance = probfail   !use probability of failure for site investigation performance
		
				else !cost metric should be 1 or greater here
		
					si_performance = diffgeo
		
				end if
 
end subroutine

!CDF of the lognormal distribution
elemental real function logcdf(x,mu,std)

	real, intent(in) ::  x,mu,std !x, mean and standard deviation of the log of values

	logcdf = 0.5+0.5*erf((log(x)-mu)/(std*1.41421))

end function

end module
