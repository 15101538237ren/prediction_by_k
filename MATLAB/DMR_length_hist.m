function DMR_length_hist()
    dmr_length_fp="../DATA/Genomic_Features/DMR/region_length_of_merged_dmr.txt";
    fig_path = "../FIGURES/DMR_length_hist.pdf";
    M=readmatrix(dmr_length_fp);
    N_BINS=20;
    fig=figure(1);
    subplot(2,2,1)
    mean_len= mean(M);
    median_len=median(M);
    mode_len= mode(M);
    [histFreq, histXout] = hist(M, N_BINS);
    bar(histXout, histFreq/sum(histFreq), 'b');
    ylim([0,1]);
    title("Mean: "+ num2str(round(mean_len, 2)) +", Median:" + num2str(round(median_len, 2)) +", Mode:" + num2str(round(mode_len, 2)));
    xlabel("length of DMR");
    ylabel("Proportion");
    print(fig, fig_path, '-dpdf','-opengl','-r300');
    close(fig);
end