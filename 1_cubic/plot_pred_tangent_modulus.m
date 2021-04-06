%% plot predicted tangent modulus
figure
% Sparse_id model
y1 = reshape(fitting_tmd(1,1,stepf),1,[]);
y2 = reshape(fitting_tmd(2,2,stepf),1,[]);
y3 = reshape(fitting_tmd(3,3,stepf),1,[]);
plot(step_range, y1,'r-','Linewidth',2); hold on;
plot(step_range, y2,'b-','Linewidth',2);
plot(step_range, y3,'g-','Linewidth',2);
if predict_path == 4
    y4 = reshape(fitting_tmd(6,6,stepf),1,[]);
    plot(step_range,y4,'m-','Linewidth',2);
end

% reference data
r1 = data{predict_path,4}(:,1);
r2 = data{predict_path,4}(:,7);
r3 = data{predict_path,4}(:,12);
plot(step_range,r1,'ro');
hold on;
plot(step_range,r2,'b^');
plot(step_range,r3,'gs');
if predict_path  == 4
    r4 = data{predict_path,4}(:,21);
    plot(step_range,r4,'md');
    legend({'D_{1111} of spid','D_{2222} of spid','D_{3333} of spid','D_{1212} of spid',...
        'D_{1111} of data','D_{2222} of data','D_{3333} of data','D_{1212} of data'},...
        'Location','southeast');
else
    
    
    legend({'D_{1111} of spid','D_{2222} of spid','D_{3333} of spid',...
        'D_{1111} of data','D_{2222} of data','D_{3333} of data'},...
        'Location','southeast');
end
grid on;

xlabel('strain');
ylabel('Tangent modulus  [GPa]');

t = title(title_name);
t.FontName = 'Tahoma';
t.FontWeight = 'bold';
t.FontSize = 13;
set(gca,'FontName','Tahoma','FontSize',13,'FontWeight','bold');
set(gcf,'color','white');
set(gcf, 'OuterPosition',[20,500,350,400])