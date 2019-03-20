function [] = manual_roi_sorting(cn)

N = cn.n_cells;
good_purkinje = zeros(N,1);
cell = 1;

fig = figure('Name','Masks');
set(fig,'KeyPressFcn',@keypress);


    function redraw() 
                clf
                hold on
                
                if cell < N
                    msk = cn.mask{1,cell};
                    imshow(msk)
                    
                    title(['ROI n° ', num2str(cell)])
                else
                    disp('All cells done.')
                end
                
     end
 
redraw()

    function keypress(~,evnt) %keypress reaction
                switch lower(evnt.Key)
                    case 'rightarrow'
                        if cell<=N
                            good_purkinje(cell) = 1; %we keep it
                            cell = cell + 1;
                        end
                    case 'leftarrow'
                        if cell<=N
                            good_purkinje(cell) = 0;
                            cell = cell + 1;
                        end
                    otherwise
                        return  %breaks it 
                end
                redraw()
           end

    global good_purkinje

end



