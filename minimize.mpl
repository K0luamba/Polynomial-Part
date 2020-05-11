#eqVec - список из правой части, полиномов по порядку степеней, comments - вывод комментариев по ходу решения
local minimize := proc(eqVec, comments := true)
	description "Минимизация правой части уравнения и выделение полиномиального слагаемого решения";
    local prevEq, leftPart, currEq, currEqVec, i, j, k, w, maxV, m, transLeft, a, subsY, needStop, polPart;
    currEqVec := eqVec;
    polPart := 0;
    currEq := 0;

    for j while currEqVec[1] <> 0 do

        # получаем из currEqVec уравнение для вывода текущего результата
        leftPart := currEqVec[2]*y(x);
        for i from 3 to nops(currEqVec) do
        	leftPart := leftPart + currEqVec[i]*diff(y(x), x$(i-2));
        end do;
        prevEq := currEq;
        currEq := leftPart = currEqVec[1];
        if comments then 
        	print("Текущее уравнение:");
        	print(currEq);
        end if;

		#проверка случая, когда в правой части есть член со степенью, меньшей 0
    	if ldegree(currEqVec[1]) < 0 or degree(currEqVec[1]) < 0 then
    		if comments then 
    			print("Алгоритм прекратил работу по 2ому условию остановки.") 
    		end if;
    		polPart := polPart - subsY; #учитываем, что последнее слагаемое полиномаильной части лишнее
    		if polPart <> 0 and comments then 
	        	print("Выделили полиномиальное слагаемое решения.");
	        end if;
    		return [subs(y=z, prevEq), polPart];
    	end if;

        # получаем величину w(омега) и m
        w := degree(currEqVec[2]); #для члена при y(x)
        for i from 3 to nops(currEqVec) do
        	maxV := degree(currEqVec[i]) - (i-2);
        	if maxV > w then
        		w := maxV;
        	end if;
        end do;
        m := degree(currEqVec[1]) - w; # m = n - w

        # проверяем случай, когда для всех указанных ниже i_j : i_j > m >= 0
        if m >= 0 then
        	needStop := true;
	        for i from 2 to nops(currEqVec) do
	        	if (degree(currEqVec[i]) - (i-2)) = w and i-2 <= m then
	        		needStop := false;
	        	end if;
	        end do;
	        if needStop then
	        	if comments then
        			print("Алгоритм прекратил работу по 1ому условию остановки.");
        		end if;
        		if polPart <> 0 and comments then
		        	print("Выделили полиномиальное слагаемое решения.");
		        end if;
        		return [subs(y=z, currEq), polPart];
        	end if;
	    end if;

        # поиск величины a
        transLeft := PDEtools:-dsubs(y(x) = x^m, leftPart);
        a := lcoeff(currEqVec[1]) / lcoeff(transLeft);

        # производим замену (т.е. изменение вектора уравнения currEqVec)
        subsY := a*x^m;
        if comments then
	        print("Производим замену:");
	        print(y(x) = subsY + z(x));
	    end if;
        polPart := polPart + subsY; #учитываем замену для составления полиномиальной части решения
	        currEqVec[1] := currEqVec[1] - subsY*currEqVec[2];
        for i from 3 to nops(currEqVec) do
        	currEqVec[1] := currEqVec[1] - diff(subsY, x$(i-2))*currEqVec[i];
        end do;
        currEqVec[1] := simplify(currEqVec[1]);
    end do;

    # получаем уравнение для вывода окончательного результата
    leftPart := currEqVec[2]*y(x);
    for i from 3 to nops(currEqVec) do
    	leftPart := leftPart + currEqVec[i]*diff(y(x), x$(i-2));
    end do;
    currEq := leftPart = currEqVec[1];
    if comments then
    	print("Получили полиномиальное решение.");
	end if;
    return [subs(y=z, currEq), polPart];
end proc;