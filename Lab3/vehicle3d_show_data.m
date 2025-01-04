function vehicle3d_show_data(t,x,u,printFigs,filename,ifig)

    % test input variables
    if ~exist('printFigs','var') | isempty(printFigs), printFigs = 0; end
    if ~exist('filename','var') | isempty(filename), filename = 'test'; end
    if ~exist('ifig','var') | isempty(ifig), ifig = 0; end

    % prepare data for plots
    if ~iscell(x)
        t = {t};
        x = {x};
        u = {u};
    end
    nD = length(x);
    nx = size(x{1},1);
    nu = size(u{1},1);
    NSim = size(x{1},2);
    has_m_dyn = 0;

    for iD = 1:nD
        p{iD} = x{iD}(1:3,:);
        v{iD} = x{iD}(4:6,:);
        if nx == 18
            for k = 1:NSim
                lbd{iD}(:,k) = R2Euler(reshape(x{iD}(7:15,k),3,3));
            end
            om{iD} = x{iD}(16:18,:);
        elseif nx == 19
            for k = 1:NSim
                lbd{iD}(:,k) = R2Euler(reshape(x{iD}(7:15,k),3,3));
            end
            om{iD} = x{iD}(16:18,:);
            m{iD} = x{iD}(19,:);
            has_m_dyn = 1;
        elseif nx == 12
            lbd{iD} = x{iD}(7:9,:);
            om{iD} = x{iD}(10:12,:);
        elseif nx == 13
            lbd{iD} = x{iD}(7:9,:);
            om{iD} = x{iD}(10:12,:);
            m{iD} = x{iD}(13,:);
            has_m_dyn = 1;
        end
    end
    
    % 3D plot
    ifig = ifig+1;
    figure(ifig);
    leg_str = '';
    for iD = 1:nD
        hplt(iD) = plot3(p{iD}(1,:),p{iD}(2,:),p{iD}(3,:));
        if iD == 1, hold on; end
        hini = plot3(p{iD}(1,1),p{iD}(2,1),p{iD}(3,1),'og','MarkerSize',2);
        hend = plot3(p{iD}(1,end),p{iD}(2,end),p{iD}(3,end),'xr','MarkerSize',2);
        leg_str = [leg_str ',''Sim' num2str(iD) ''''];
    end
    hold off;
    grid on;
    axis equal;
    xlabel('$$p_1$$ [m]');
    ylabel('$$p_2$$ [m]');
    zlabel('$$p_3$$ [m]');
    set(gca,'Zdir','reverse');
    set(gca,'Ydir','reverse');
    title('Vehicle trajectory');
    eval(['legend([hplt,hini,hend]' leg_str ',''ini'',''end'');']);
    print2pdf([filename '_traj'],printFigs);
    
    % actuation
    ifig = ifig+1;
    figure(ifig);
    for iu = 1:nu
        subplot(nu,1,iu);
        for iD = 1:nD
            plot(t{iD},u{iD}(iu,:));
            if iD == 1, hold on; end
        end
        hold off;
        grid on;
        ylabel(['$$u_' num2str(iu) '(t)$$']);
    end
    xlabel('$$t$$ [s]');
    title('Control variables');
    print2pdf([filename '_act'],printFigs);

    ifig = ifig+1;
    figure(ifig);
    for ivar = 1:3
        subplot(3,1,ivar);
        for iD = 1:nD
            plot(t{iD},p{iD}(ivar,:));
            if iD == 1, hold on; end
        end
        hold off;
        grid on;
        ylabel(['$$p_' num2str(ivar) '(t) [m]$$']);
    end
    xlabel('$$t$$ [s]');
    title('Vehicle position');
    print2pdf([filename '_pos'],printFigs);

    ifig = ifig+1;
    figure(ifig);
    subplot(311);
    for ivar = 1:3
        subplot(3,1,ivar);
        for iD = 1:nD
            plot(t{iD},v{iD}(ivar,:));
            if iD == 1, hold on; end
        end
        hold off;
        grid on;
        ylabel(['$$v_' num2str(ivar) '(t) [m/s]$$']);
    end
    xlabel('$$t$$ [s]');
    title('Linear velocity');
    print2pdf([filename '_vel'],printFigs);
    
    ifig = ifig+1;
    figure(ifig);
    subplot(311);
    for ivar = 1:3
        subplot(3,1,ivar);
        for iD = 1:nD
            plot(t{iD},lbd{iD}(ivar,:)*180/pi);
            if iD == 1, hold on; end
        end
        hold off;
        grid on;
        ylabel(['$$\lambda_' num2str(ivar) '(t) [deg]$$']);
    end
    xlabel('$$t$$ [s]');
    title('Attitude (euler angles)');
    print2pdf([filename '_eul'],printFigs);
    
    ifig = ifig+1;
    figure(ifig);
    for ivar = 1:3
        subplot(3,1,ivar);
        for iD = 1:nD
            plot(t{iD},om{iD}(ivar,:)*180/pi);
            if iD == 1, hold on; end
        end
        hold off;
        grid on;
        ylabel(['$$\omega_' num2str(ivar) '(t) [deg/s]$$']);
    end
    xlabel('$$t$$ [s]');
    title('Angular velocity');
    print2pdf([filename '_om'],printFigs);

    if has_m_dyn
        ifig = ifig+1;
        figure(ifig);
        for iD = 1:nD
            plot(t{iD},m{iD}(1,:));
            if iD == 1, hold on; end
        end
        grid on;
        xlabel('$$t$$ [s]');
        ylabel('$$m(t)$$ [kg]');
        title('Mass');
        print2pdf([filename '_mass'],printFigs);
    end

end