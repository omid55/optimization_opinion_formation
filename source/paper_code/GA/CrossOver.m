%Omid55
function [ population,fitnesses ] = CrossOver( population,fitnesses,Pc,o )

%disp('CrossOver ... ');
%% Calculation Body
% parents = SelectParents(fitnesses);
% for i=1:length(parents)/2
%     
%     p1 = parents(i);
%     find(parents ~= p1);
%     
% end

len = size(population,1);
children = [];
for i=1:len
    
    fit = fitnesses;
    p1 = RouletteSelect(fit);
    fit(p1) = [];
    p2 = RouletteSelect(fit);
    if p2>p1
        p2 = p2 + 1;
    end
    
    prob = rand(1,1);
    if prob <= Pc
        % do crossover
        
        % 3 part crossover
        list = 1:size(population,2);
        i1 = randi(length(list),1);
        list(i1) = [];
        if i1+1 <= length(list)
            list(i1+1) = [];
        end
        if i1-1 >= 1
            list(i1-1) = [];
        end
        i2 = randi(length(list),1);
        i2 = list(i2);
        mx = max(i1,i2);
        mi = min(i1,i2);
        i1 = mi;
        i2 = mx;
        
        part = randi(3,1);
        
        if part == 1
            children = [children; doCrossover(population(p2,:),population(p1,:),i1)];
            children = [children; doCrossover(population(p1,:),population(p2,:),i1)];
        else
            if part == 2
                children = [children; doCrossover(population(p1,:),population(p2,:),i1,i2)];
                children = [children; doCrossover(population(p2,:),population(p1,:),i1,i2)];
            else  % it means part == 3
                children = [children; doCrossover(population(p2,:),population(p1,:),i2)];
                children = [children; doCrossover(population(p1,:),population(p2,:),i2)];
            end
        end
        
    end
    
end

childrenFitnesses = CalculateFitnesses(children,o);
population = [population; children];
fitnesses = [fitnesses; childrenFitnesses];


end
