subroutine priorvalue(point,water_depth,pv_min,pv_max)

	real point,water_depth,pv_min,pv_max,v

        if(point <= water_depth) then
                pv_min=-5.0
                pv_max=-2.0

        else if(point <=13) then
                pv_min=2.0
                pv_max=4.0 

        else if(point <=17.) then
                pv_min=3.2 
                pv_max=4.5

        else if(point <=30.) then
                pv_min=3.5
                pv_max=4.8 

        else
                pv_min=4.1 
                pv_max=4.9 
        endif

return
end
