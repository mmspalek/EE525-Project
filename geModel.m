\section*{Appendix D: getModel.m}

\begin{lstlisting}[style=Matlab-editor]

function [a_model, p_model] = getModel(T, horizon)
    t = 0:T:horizon;
    a_model = zeros(1, length(t));
    p_model = zeros(1, length(t));
    poi = [ % acceleration up
            02.60 03.91 00.000 00.242;
            03.91 04.69 00.242 00.242;
            04.69 06.00 00.242 00.000;
            % de-acceleration up
            09.69 10.90 00.000 -0.425;
            10.90 12.07 -0.425 00.000;
            % acceleration dn
            24.00 25.30 00.000 -0.337;
            25.30 25.80 -0.337 -0.337;
            25.80 26.60 -0.337 00.000;
            % de-acceleration dn
            30.65 31.70 00.00 00.445;
            31.70 33.00 00.445 00.00];
    for i = 1:1:height(poi)
        fro = poi(i, 1); to = poi(i, 2);
        init = poi(i, 3); fin = poi(i, 4);
        range = round(fro/T)+1:round(to/T)+1;
        a_model(range) = (((t(range) - fro) / (to - fro)) * (fin - init)) + init;
    end
    for i = 1:1:length(t)
        if (t(i) >= 12) && (t(i) <= 24)
            p_model(i) = 3.33;
        end
    end
end

\end{lstlisting}