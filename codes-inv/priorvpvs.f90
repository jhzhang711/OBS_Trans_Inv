subroutine priorvpvs(point,v,water_depth,pvpvs_min,pvpvs_max)

	real point,v,water_depth
	real pvpvs_min,pvpvs_max


        if(point <= water_depth) then
		pvpvs_min=-999
		pvpvs_max=-999

	else if (point<=13.0) then
		pvpvs_min=1.55 
		pvpvs_max=2.20

	else if (point<=45.0) then
		pvpvs_min=1.65 
		pvpvs_max=2.20 

!jiahui try

        else
		pvpvs_min=1.60 
		pvpvs_max=2.20
        endif


return
end
