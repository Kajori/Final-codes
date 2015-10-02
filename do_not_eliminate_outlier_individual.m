%here the malacious node sends the same value to all the other nodes
%Remove min-max of each eigen value and replace with the average of the
%values
function [neighbour_val] = eliminate_outlier_individual(x,j,NO_AREA,deg,neighbour)
%neighbour = index of the neighbour
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %         find the top f and bottom f outliers
    %         max_sum and min_sum are the top and bottom outliers
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
    index_i_values=zeros(1,NO_AREA); %NO_AREA because one does not use one's own row
    neighbour_val=zeros(deg,1);
       
    for i=1:deg
        for area=1:NO_AREA
                index_i_values(1,area)=x((area-1)*deg+i,j);
        end %end of col=1:NO_AREA
        [max_val, max_index ]=max(index_i_values);
        [min_val, min_index ]=min(index_i_values);
           
        replace=0;
        if(neighbour==max_index)%Find replaces
            replace=1;
            for col=1:NO_AREA
                if(max_val==x((col-1)*deg+i,j))
                    replace=0;
                end
            end
        elseif(neighbour==min_index)  
            replace=1;
            for col=1:NO_AREA 
                 if(min_val==x((col-1)*deg+i,j))
                     replace=0;
                 end
            end
        end
        if(replace==0)
             neighbour_val(i)=x((neighbour-1)*deg+i,j);
        else
             neighbour_val(i)=mean(index_i_values);
        end
    end %end for i=1 to deg
end
        
