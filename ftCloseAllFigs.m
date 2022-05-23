function ftCloseAllFigs()

        fprintf('Closing open figures...\n');
        figHandles = get(groot, 'Children');
        if (~isempty(figHandles))
            close(figHandles);
        end