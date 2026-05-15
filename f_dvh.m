function R = f_dvh(dose_volume, roi_masks, roi_names, roi_colors, voxel_vol, roi_tumor)
    n_roi = length(roi_masks);
    R = struct('name', cell(n_roi,1), 'color', cell(n_roi,1), ...
               'dvh_dose', cell(n_roi,1), 'dvh_vol', cell(n_roi,1), ...
               'vol_cc', cell(n_roi,1), ...
               'Dmax', cell(n_roi,1), 'Dmean', cell(n_roi,1), 'Dmin', cell(n_roi,1), ...
               'D98', cell(n_roi,1), 'D95', cell(n_roi,1), 'D50', cell(n_roi,1), 'D2', cell(n_roi,1), ...
               'V20', cell(n_roi,1), 'V10', cell(n_roi,1), ...
               'tumor', cell(n_roi,1));

    for i = 1:n_roi
        mask = roi_masks{i};
        dvh_dose = dose_volume(mask);
        dvh_dose = dvh_dose(~isnan(dvh_dose));
        dvh_dose = dvh_dose(:);

        if isempty(dvh_dose)
            R(i).name = roi_names{i};
            R(i).color = roi_colors(i,:);
            R(i).vol_cc = 0;
            R(i).dvh_dose = [];
            R(i).dvh_vol = [];
            R(i).Dmax = NaN; R(i).Dmean = NaN; R(i).Dmin = NaN;
            R(i).D98 = NaN; R(i).D95 = NaN; R(i).D50 = NaN; R(i).D2 = NaN;
            R(i).V20 = NaN; R(i).V10 = NaN;
            R(i).tumor = roi_tumor(i);
            continue;
        end

        vol_cc = voxel_vol * length(dvh_dose);

        dvh_dose_sorted = sort(dvh_dose, 'ascend');
        dvh_vol = 100 * (1 - (0:length(dvh_dose_sorted)-1)' / length(dvh_dose_sorted));

        Dmax = max(dvh_dose);
        Dmean = mean(dvh_dose);
        Dmin = min(dvh_dose);

        if length(dvh_dose_sorted) >= 3
            D98 = interp1(dvh_vol, dvh_dose_sorted, 98, 'linear', 'extrap');
            D95 = interp1(dvh_vol, dvh_dose_sorted, 95, 'linear', 'extrap');
            D50 = interp1(dvh_vol, dvh_dose_sorted, 50, 'linear', 'extrap');
            D2  = interp1(dvh_vol, dvh_dose_sorted, 2,  'linear', 'extrap');
        else
            D98 = Dmin; D95 = Dmin; D50 = Dmean; D2 = Dmax;
        end

        V20 = interp1(dvh_dose_sorted, dvh_vol, 20, 'linear', 'extrap');
        V10 = interp1(dvh_dose_sorted, dvh_vol, 10, 'linear', 'extrap');

        R(i).name = roi_names{i};
        R(i).color = roi_colors(i,:);
        R(i).vol_cc = vol_cc;
        R(i).dvh_dose = dvh_dose_sorted;
        R(i).dvh_vol = dvh_vol;
        R(i).Dmax = Dmax; R(i).Dmean = Dmean; R(i).Dmin = Dmin;
        R(i).D98 = D98; R(i).D95 = D95; R(i).D50 = D50; R(i).D2 = D2;
        R(i).V20 = V20; R(i).V10 = V10;
        R(i).tumor = roi_tumor(i);
    end
end
