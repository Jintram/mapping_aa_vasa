
give_better_textsize_plot_tmp = 
    function(TEXTSIZE, myFamily='Arial'){
      
      if (is.null(myFamily)) {
        theme(#legend.position="none",
              text = element_text(size=TEXTSIZE),
              axis.text = element_text(size=TEXTSIZE),
              plot.title = element_text(size=TEXTSIZE))
      } else {
        print('Set to Arial')
        theme(#legend.position="none",
              text = element_text(size=TEXTSIZE, family=myFamily),
              axis.text = element_text(size=TEXTSIZE, family=myFamily),
              plot.title = element_text(size=TEXTSIZE, family=myFamily))
      }
    
    }




samples_read_counts = c(4000, 3020520,  6041040,  9061560, 12082080, 15102600, 18123120, 21143640, 24164160, 27184680, 30205200)

dir_path = '/Users/m.wehrens/Data/2022_09_micetimeline/2022_11_mice37samples/countables_ds-test/'

collected_sums =
    sapply(samples_read_counts, function(cnt) {
        
            data_table = read.table(paste0(dir_path, 'head',cnt,'__Martijn-Wehrens-sample1_AACGFYKM5_S2_L001_pT.nonRibo_E99_Aligned.out.counts.tsv'), header=1)
            # View(data_table)
            
            return(sum(data_table$count))
            
        })

df_results_ds_mapping =
    data.frame(count_sum=collected_sums, input_lines=samples_read_counts, input_reads=samples_read_counts/4)


lin_fit_end = lm(df_results_ds_mapping$count_sum[6:11]/1e6 ~ poly(df_results_ds_mapping$input_reads[6:11]/1e6, 1, raw = T))
lin_fit_end_coeffs=lin_fit_end$coefficients

p=ggplot(df_results_ds_mapping, aes(x=input_reads/1e6, y=count_sum/1e6))+
    geom_point()+
    geom_line()+theme_bw()+
    geom_hline(yintercept = max(df_results_ds_mapping$count_sum)/1e6)+
    geom_vline(xintercept = 7551300/1e6)+
    geom_hline(yintercept = 0)+
    geom_vline(xintercept = 0)+
    geom_hline(yintercept = max(df_results_ds_mapping$count_sum)/1e6*2)+
    geom_vline(xintercept = 7551300/1e6*2)+
    geom_abline(slope = lin_fit_end_coeffs[2], intercept = lin_fit_end_coeffs[1], linetype = "dashed")+
    geom_abline(slope = 1, intercept = 0, linetype = "dashed")+xlim(c(0,7551300*2/1e6))+ylim(c(0,max(df_results_ds_mapping$count_sum)/1e6*2))+
    give_better_textsize_plot_tmp(9)+
    ggtitle(paste0('fit based last points; y0=', round(lin_fit_end_coeffs[1],3), '; a=', round(lin_fit_end_coeffs[2],3),'\n',
                   'max(x)=',round(max(df_results_ds_mapping$input_reads)/1e6,2),'; max(y)=', round(max(df_results_ds_mapping$count_sum)/1e6,3)))

p
ggsave(plot=p, filename = '/Users/m.wehrens/Data/2022_09_micetimeline/2022_11_mice37samples/countables_ds-test/plot.pdf',
       device = cairo_pdf, width = 75, height = 75,units = 'mm')






