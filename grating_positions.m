% Calculate grating positions to cover whole range from 4000 A to 8000 A with given overlap
overlap = 20;
range1st = [5100, 8000];
range2nd = [4000, 5100];

% First Order: [Increment, Degrees, Start_A, End_A]
raw1st = [
	5321, 35.78, 8401.4, 8869.9;
	5430, 35.25, 8197.2, 8666.3;
	6382, 30.62, 6406.2, 6878.6;
	6457, 30.25, 6263.5, 6736.0;
	6943, 27.88, 5328.6, 5801.4
	];

% Second Order: [Increment, Degrees, Start_A, End_A]
raw2nd = [
	4712, 38.75, 4758.8, 4991.1;
	5040, 37.15, 4459.9, 4693.4;
	5242, 36.17, 4273.7, 4507.7
	];


process_order('FIRST ORDER', raw1st, range1st(1), range1st(2), overlap);
process_order('SECOND ORDER', raw2nd, range2nd(1), range2nd(2), overlap);

function process_order(name, data, orderStart, orderEnd, overlap)
fprintf('\n========================================================================\n');
fprintf(' CALCULATION FOR: %s  (Target: %.0f - %.0f A)\n', name, orderStart, orderEnd);
fprintf('========================================================================\n');

incVals   = data(:, 1);
degVals   = data(:, 2);
startVals = data(:, 3);
endVals   = data(:, 4);

centers = (startVals + endVals) / 2;

widths = endVals - startVals;
avgWidth = mean(widths);

modelInc = polyfit(centers, incVals, 1);

% We use the raw calibration data to map motor steps to physical degrees
% TODO: there surely must be the formula listed somewhere, or I could just
% solve the linear equation system myself to get it, but it doesn't matter 
% for now
modelAng = polyfit(incVals, degVals, 1);



fprintf('Avg Window Size: %.0f A\n', avgWidth);
fprintf('Model: Increment = %.4f * WL + %.4f\n', modelInc(1), modelInc(2));
fprintf('Model: Angle     = %.4f * Inc + %.4f\n', modelAng(1), modelAng(2));
fprintf('------------------------------------------------------------------------\n');
fprintf('%-10s | %-10s | %-15s | %-20s\n', 'Increment', 'Angle(deg)', 'Center WL (A)', 'Approx Range (A)');
fprintf('------------------------------------------------------------------------\n');

currentCenter = orderStart + (avgWidth / 2);

while (currentCenter - avgWidth/2) < orderEnd
	calcInc = polyval(modelInc, currentCenter);
	calcAng = polyval(modelAng, calcInc);
	currentStart = currentCenter - (avgWidth/2);
	currentEnd   = currentCenter + (avgWidth/2);
	fprintf('%-10.0f | %-10.2f | %-15.1f | %.1f - %.1f\n', ...
		calcInc, calcAng, currentCenter, currentStart, currentEnd);

	currentCenter = currentCenter + (avgWidth - overlap);
	if currentCenter > (orderEnd + 2000)
		break;
	end
end
end
