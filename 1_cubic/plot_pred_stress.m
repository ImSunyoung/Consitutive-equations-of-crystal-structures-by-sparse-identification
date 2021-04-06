%% plot predicted second PK stress 
title_name = sprintf('Cubic : %s', materials{type});

[~, tip] = max(abs(data{predict_path,2}(end,:)));
step_range = data{predict_path,2}(:,tip);
stepf = [1:nStep];


% Sparse_id model
y1 = reshape(fitting_PK2(1,1,stepf),1,[]);
y2 = reshape(fitting_PK2(2,2,stepf),1,[]);
y3 = reshape(fitting_PK2(3,3,stepf),1,[]);
plot(step_range, y1,'r:','Linewidth',2);
hold on;
plot(step_range, y2,'b:','Linewidth',2);
plot(step_range, y3,'g:','Linewidth',2);
if predict_path == 4
    y4 = reshape(fitting_PK2(1,2,stepf),1,[]);
    plot(step_range,y4,'m:','Linewidth',2);
end

% reference data
r1 = data{predict_path,3}(:,1);
r2 = data{predict_path,3}(:,2);
r3 = data{predict_path,3}(:,3);
plot(step_range,r1,'ro');
hold on;
plot(step_range,r2,'b^');
plot(step_range,r3,'gs');
if predict_path == 4
    r4 = data{predict_path,3}(:,6);
    plot(step_range,r4,'md');
    legend({'PK2 11_{spid}','PK2 22_{spid}','PK2 33_{spid}','PK2 12_{spid}',...
        'PK2 11_{data}','PK2 22_{data}','PK2 33_{data}','PK2 12_{data}'},...
        'Location','southeast');
else
    legend({'PK2 11_{spid}','PK2 22_{spid}','PK2 33_{spid}',...
        'PK2 11_{data}','PK2 22_{data}','PK2 33_{data}'},...
        'Location','southeast');
end
grid on;
xlabel('strain');
ylabel('stress   [GPa]');
t = title(title_name);
t.FontName = 'Tahoma';
t.FontWeight = 'bold';
t.FontSize = 13;
set(gca,'FontName','Tahoma','FontSize',13,'FontWeight','bold');
set(gcf,'color','white');
set(gcf, 'OuterPosition',[20,500,350,400])