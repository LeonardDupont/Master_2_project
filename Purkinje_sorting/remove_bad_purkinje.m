function [cn_bis] = remove_bad_purkinje(cn,good_purkinje)

  
    N = cn.n_cells;
    disp(['Removing fake rois based on input vector. Initial number of cells = ',num2str(N)])
    cn_bis.n_cells = sum(good_purkinje);
    cn_bis.fov_height = cn.fov_height;
    cn_bis.fov_width = cn.fov_width;
    
    cell = 1;
    
    for k = 1:N
        if good_purkinje(k)
            % . . . . create fields again . . . . . . . . . . . . . . . . .
            cn_bis.intensity(:,cell) = cn.intensity(:,k);
            cn_bis.mask{1,cell} = cn.mask{1,k};
            cn_bis.roi{1,cell} = cn.roi{1,k};
            cn_bis.roi_landscape{1,cell} = cn.roi_landscape{1,k};
            cn_bis.centroid{1,cell} = cn.centroid{1,k};
            if isfield(cn,'spikes')
                cn_bis.spikes(:,cell) = cn.spikes(:,k);
            end
            % . . . . . . . . . . . . . . . . .. . . . . . . . . . . . . . 
            cell = cell + 1;
        end
    end
    
    disp(['Done, final number of cells = ',num2str(cn_bis.n_cells),'.'])
end