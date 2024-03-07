

library(tidyverse)
library(mgcv)
library(metafor)
 
df <- read_csv("smid_data.csv") %>%
  mutate(n_subjects = n_exp + n_control) %>%
  mutate(total_iron_taken_mg = elemental_iron_dose_mg_per_day*intervention_length_weeks*7)


df
df %>% glimpse()


df %>% 
  ggplot(aes(x=initial_ferritin_ng_ml , y=ferritin_effect_size_smd, size=n_subjects)) + 
  geom_point() + 
  lims(x=c(0,50), y=c(-2,6))


# --- use metafor to reproduce meta-analysis


meta <- rma(ferritin_effect_size_smd, sei = ferritin_std_error_smd, slab=id, 
            data = df, method="REML", measure="SMD")
meta
#Results agree, modulo minor differences in 95% CI calculations


# --- Make forest plot ----
forest(meta, header="Study", mlab="Overall pooled results", cex = 1.5)

# Generate the forest plot PNG
png("Meta-analysis forest plot of iron supplements and ferritin level increase.png", width=1200, height=600)
forest(meta, header="Study", mlab="Overall pooled results", cex=1.75,
       main="Effect of iron supplements on ferritin level",
       cex.main = 2.5)

dev.off()



#Plotting setup 
xs <- seq(5, 50, length=500)


# --- Linear fit, purely for plotting purposes --- 



res.lin <- rma(ferritin_effect_size_smd, sei = ferritin_std_error_smd, slab=id, 
               mods = ~ initial_ferritin_ng_ml,
               data = df, method="REML", measure="SMD")
res.lin


#Linear plot
regplot(res.lin, las=1, digits=1, bty="l", psize=.75/df$ferritin_std_error_smd,
        xlab="Predictor", main="Linear Model")


# --- MGCV fit ---

K <- 4

library(mgcv)
sm <- smoothCon(s(initial_ferritin_ng_ml, bs="cr", k=K), data=df, absorb.cons=TRUE)[[1]]
res.tps <- rma(ferritin_effect_size_smd, sei = ferritin_std_error_smd, slab=id, 
               mods = ~ sm$X,
               data = df, method="REML", measure="SMD")


res.tps


sav <- predict(res.tps, newmods=PredictMat(sm, data.frame(initial_ferritin_ng_ml=xs)))
regplot(res.lin, mod=2, pred=sav, xvals=xs, las=1, digits=1, bty="l",
        psize=.75/df$ferritin_std_error_smd, 
        xlab="Predictor", main="Thin Plate Spline Model")


met_plot <- data.frame(ferritin = xs,
                      yhat = sav$pred,
                      ci_lo = sav$ci.lb,
                      ci_hi = sav$ci.ub)



# -- Setup plot --


point_fill <- "#deebf7"
point_color <- "black"

lwd <- 0.5
draw_lwd <- 0.15
fnt <- 8
title_fnt <- 8
subtitle_fnt <- 6
wm_fnt <- 7
hline_lwd <- 0.25
grid_lwd <- 0.2

alf <- 0.12


bubble_max <- 3.5

#, range=c(0.5,4)

plt <- met_plot %>%
  ggplot(aes(x=ferritin, y=yhat)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = hline_lwd) + 
  geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi), fill = "#3182bd", alpha = alf) + 
  geom_line(color = "#3182bd", linewidth = lwd) + 
  #Study bubbles
  geom_point(aes(x=initial_ferritin_ng_ml , y=ferritin_effect_size_smd, size=2/ferritin_std_error_smd ),
             data = df, 
             shape=21,
             fill = point_fill,
             stroke = 0.2,
             color = point_color) + 
  
  #Change point sizing
  scale_size(guide="none", range=c(0.5,bubble_max)) +
  #Set axis limits
  scale_x_continuous(limits = c(0,50), expand = c(0,0),
                     breaks = seq(0,50, by=5),
                     name = "Ferritin level (ng/mL)") + 
  scale_y_continuous(limits = c(-1,8), expand = c(0,0),
                     breaks = seq(-2,8,by=1),
                     name = "Ferritin increase (effect size)") + 
  coord_cartesian(ylim = c(-1,6.5))+
  ggtitle(label = "Iron supplements have dramatic effects when ferritin is <20 ng/mL", subtitle = "Data from Smid et al. 2024") + 
  labs(caption = "RunningWritings.com") + 
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5, size=title_fnt),
        plot.caption.position = "panel",
        plot.caption = element_text(face = "bold", size = wm_fnt, vjust=0),
        plot.subtitle = element_text(hjust = 0.5, size=subtitle_fnt),
        panel.grid.major = element_line(linewidth = grid_lwd),
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = fnt, color = "black"),
        axis.text = element_text(size = fnt, color = "black"))

plt


ggsave("Iron supplement and ferritin level increase - effect sizes.png", plot=plt,
       width = 1200, height = 800, units = "px")





# --- Simple plot for top of blog post


lwd <- 0.75


plt <- met_plot %>%
  ggplot(aes(x=ferritin, y=yhat)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = hline_lwd) + 
  geom_line(color = "#3182bd", linewidth = lwd) + 
  #Change point sizing
  scale_size(guide="none", range=c(0.5,bubble_max)) +
  #Set axis limits
  scale_x_continuous(limits = c(0,50), expand = c(0,0),
                     breaks = seq(0,50, by=5),
                     name = "Ferritin level (ng/mL)") + 
  scale_y_continuous(limits = c(-1,8), expand = c(0,0),
                     breaks = seq(-2,8,by=1),
                     name = "Ferritin increase (effect size)") + 
  coord_cartesian(ylim = c(-1,6.5))+
  ggtitle(label = "Iron supplements have dramatic effects when ferritin is <20 ng/mL") + 
  labs(caption = "RunningWritings.com") + 
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5, size=title_fnt),
        plot.caption.position = "panel",
        plot.caption = element_text(face = "bold", size = wm_fnt, vjust=0),
        plot.subtitle = element_text(hjust = 0.5, size=subtitle_fnt),
        panel.grid.major = element_line(linewidth = grid_lwd),
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = fnt, color = "black"),
        axis.text = element_text(size = fnt, color = "black"))

plt


ggsave("Iron supplement and ferritin level increase - simple.png", plot=plt,
       width = 1200, height = 800, units = "px")



# --- Bubbles only


lwd <- 0.75


plt <- met_plot %>%
  ggplot(aes(x=ferritin, y=yhat)) +
  #Study bubbles
  geom_point(aes(x=initial_ferritin_ng_ml , y=ferritin_effect_size_smd, size=2/ferritin_std_error_smd ),
             data = df, 
             shape=21,
             fill = point_fill,
             stroke = 0.2,
             color = point_color) + 
  
  #Change point sizing
  scale_size(guide="none", range=c(0.5,bubble_max)) +
  #Set axis limits
  scale_x_continuous(limits = c(0,50), expand = c(0,0),
                     breaks = seq(0,50, by=5),
                     name = "Ferritin level (ng/mL)") + 
  scale_y_continuous(limits = c(-1,8), expand = c(0,0),
                     breaks = seq(-2,8,by=1),
                     name = "Ferritin increase (effect size)") + 
  coord_cartesian(ylim = c(-1,6.5))+
  ggtitle(label = "Data from 13 iron supplement studies totaling 449 athletes", subtitle = "Data from Smid et al. 2024") + 
  labs(caption = "RunningWritings.com") + 
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5, size=title_fnt),
        plot.caption.position = "panel",
        plot.caption = element_text(face = "bold", size = wm_fnt, vjust=0),
        plot.subtitle = element_text(hjust = 0.5, size=subtitle_fnt),
        panel.grid.major = element_line(linewidth = grid_lwd),
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = fnt, color = "black"),
        axis.text = element_text(size = fnt, color = "black"))

plt

ggsave("Iron supplement and ferritin level increase - bubble plot.png", plot=plt,
       width = 1200, height = 800, units = "px")




# --------------- VO2max data ----------------





# --- use metafor to reproduce meta-analysis


meta_vo2 <- rma(vo2max_effect_size, sei = vo2max_std_error, slab=id, 
            data = df, method="REML", measure="SMD")
meta_vo2
#Results agree, modulo minor differences in 95% CI calculations

plot(meta_vo2)


# --- Make forest plot ----
forest(meta_vo2, header="Study", mlab="Overall pooled results", cex = 1.5)

# Generate the forest plot PNG
png("Meta-analysis forest plot of iron supplements and VO2max increase.png", width=1200, height=600)
forest(meta_vo2, header="Study", mlab="Overall pooled results", cex=1.75, main="Effect of iron supplements on VO2max",
       cex.main = 2.5)
dev.off()



# --- linear and spline fits -----


df_vo2 <- df %>%
  drop_na(vo2max_effect_size)
#Don't use all data for knots, just studies with VO2 measurements


#Plotting setup 
xs <- seq(5, 35, length=500) #no VO2 data for the higher ranges


res.lin_vo2 <- rma(vo2max_effect_size, sei = vo2max_std_error, slab=id, 
               mods = ~ initial_ferritin_ng_ml,
               data = df_vo2, method="REML", measure="SMD")
res.lin_vo2


#Linear plot
regplot(res.lin_vo2, las=1, digits=1, bty="l", psize=.75/df_vo2$vo2max_std_error,
        xlab="Predictor", main="Linear Model")


# --- MGCV fit ---

K <- 3 #Fewer studies, fewer knots


library(mgcv)
sm <- smoothCon(s(initial_ferritin_ng_ml, bs="cr", k=K), data=df_vo2, absorb.cons=TRUE)[[1]]
res.tps_vo2 <- rma(vo2max_effect_size, sei = vo2max_std_error, slab=id, 
               mods = ~ sm$X,
               data = df_vo2, method="REML", measure="SMD")
res.tps_vo2


sav <- predict(res.tps_vo2, newmods=PredictMat(sm, data.frame(initial_ferritin_ng_ml=xs)))
regplot(res.lin_vo2, mod=2, pred=sav, xvals=xs, las=1, digits=1, bty="l",
        psize=.75/df_vo2$vo2max_std_error, 
        xlab="Predictor", main="Thin Plate Spline Model")


met_plot_vo2 <- data.frame(ferritin = xs,
                       yhat = sav$pred,
                       ci_lo = sav$ci.lb,
                       ci_hi = sav$ci.ub)




# ---- VO2max plot 



# -- Setup plot --


point_fill <- "#efedf5"
point_color <- "black"

line_col <- "#756bb1"
rib_col <- "#756bb1"

lwd <- 0.5
draw_lwd <- 0.15
fnt <- 8
title_fnt <- 8
subtitle_fnt <- 6
wm_fnt <- 7
hline_lwd <- 0.25
grid_lwd <- 0.2

alf <- 0.15


bubble_max <- 3.5

#, range=c(0.5,4)

plt <- met_plot_vo2 %>%
  ggplot(aes(x=ferritin, y=yhat)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = hline_lwd) + 
  geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi), fill = rib_col, alpha = alf) + 
  geom_line(color = line_col, linewidth = lwd) + 
  #Study bubbles
  geom_point(aes(x=initial_ferritin_ng_ml , y=vo2max_effect_size, size=2/vo2max_std_error),
             data = df_vo2, 
             shape=21,
             fill = point_fill,
             stroke = 0.2,
             color = point_color) + 
  #Change point sizing
  scale_size(guide="none", range=c(0.5,bubble_max)) +
  #Set axis limits
  scale_x_continuous(limits = c(0,35), expand = c(0,0),
                     breaks = seq(0,35, by=5),
                     name = "Ferritin level (ng/mL)") + 
  scale_y_continuous(limits = c(-2,4), expand = c(0,0),
                     breaks = seq(-2,4,by=1),
                     name = "VO2max increase (effect size)") + 
  coord_cartesian(ylim = c(-1,4))+
  ggtitle(label = "Iron supplements increase VO2max when ferritin is <15 ng/mL", subtitle = "Data from Smid et al. 2024") + 
  labs(caption = "RunningWritings.com") + 
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5, size=title_fnt),
        plot.caption.position = "panel",
        plot.caption = element_text(face = "bold", size = wm_fnt, vjust=0),
        plot.subtitle = element_text(hjust = 0.5, size=subtitle_fnt),
        panel.grid.major = element_line(linewidth = grid_lwd),
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = fnt, color = "black"),
        axis.text = element_text(size = fnt, color = "black"))

plt


ggsave("Iron supplement and VO2max increase - effect sizes.png", plot=plt,
       width = 1200, height = 800, units = "px")


