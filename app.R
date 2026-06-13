library(shiny)
library(dplyr)
library(ggplot2)
library(DT)
library(geepack)

load("app_group_summary.RData")
load("app_regression_data.RData")

adjustment_label <- paste(
  "Patient-clustered logistic GEE with an exchangeable working correlation",
  "structure and robust standard errors; adjusted for female age, female BMI,",
  "filled endometrial preparation method, filled embryo derivation, and",
  "preparation-day endometrial thickness."
)

grade_choices <- sort(unique(app_group_summary$Embryo_morphology_clean))
day_choices <- c("D5", "D6")
stage_choices <- as.character(3:6)
cell_choices <- c("A", "B", "C")

find_group <- function(day, stage, icm, te) {
  paste0(day, "_", stage, icm, te)
}
default_groups <- app_group_summary |>
  arrange(desc(n), group) |>
  distinct(group, .keep_all = TRUE) |>
  slice_head(n = 5) |>
  pull(group)
group_part <- function(i, part) {
  g <- default_groups[pmin(i, length(default_groups))]
  if (part == "day") sub("_.*", "", g) else substr(sub(".*_", "", g), part, part)
}

ui <- fluidPage(
  tags$head(
    tags$title("SingleBlastoFET"),
    tags$style(HTML("
      body { background: #f7f8f8; color: #26302f; }
      .app-shell { max-width: 1360px; margin: 0 auto; padding: 18px 20px 28px; }
      .topbar { display:flex; align-items:flex-end; justify-content:space-between; gap:16px; margin-bottom:14px; }
      .brand h1 { font-size: 24px; margin: 0; font-weight: 700; letter-spacing: 0; }
      .brand div { font-size: 13px; color: #5d6866; margin-top: 3px; }
      .panel { background:#ffffff; border:1px solid #dfe5e3; border-radius:8px; padding:14px; }
      .control-panel { min-height: 640px; }
      .embryo-row { border-top:1px solid #edf0ef; padding-top:10px; margin-top:10px; }
      .embryo-row:first-child { border-top:0; padding-top:0; margin-top:0; }
      .btn-primary { background:#2f7d6d; border-color:#2f7d6d; }
      .metric-grid { display:grid; grid-template-columns: repeat(4, minmax(120px,1fr)); gap:10px; margin-bottom:12px; }
      .metric { background:#ffffff; border:1px solid #dfe5e3; border-radius:8px; padding:12px; }
      .metric .label { color:#667370; font-size:12px; }
      .metric .value { color:#1d2c2a; font-size:22px; font-weight:700; margin-top:2px; }
      .recommend { background:#eef3f0; border:1px solid #8fa39a; padding:14px; margin-bottom:12px; }
      .recommend .value { font-size:18px; font-weight:700; color:#365f55; }
      .note { color:#667370; font-size:12px; line-height:1.45; }
      .rank-badge { font-weight:700; color:#2f7d6d; }
      .tab-content { padding-top: 12px; }
      @media (max-width: 900px) {
        .metric-grid { grid-template-columns: repeat(2, minmax(120px,1fr)); }
      }
    "))
  ),
  div(class = "app-shell",
      div(class = "topbar",
          div(class = "brand",
              h1("SingleBlastoFET"),
              div("Non-PGT single blastocyst FET decision support based on center-specific live birth outcomes")
          )
      ),
      fluidRow(
        column(
          width = 3,
          div(class = "panel control-panel",
              selectInput("candidate_n", "Candidate embryos", choices = 2:5, selected = 2),
              uiOutput("candidate_controls"),
              actionButton("compare", "Compare", class = "btn btn-primary", width = "100%"),
              tags$hr(),
              div(class = "note",
                  "Ranking prioritizes adjusted predicted live birth probability. Displayed observed rates remain center-specific descriptive estimates; sparse groups require caution.")
          )
        ),
        column(
          width = 9,
          div(class = "recommend", div(class = "label", "Recommended candidate embryo"), div(class = "value", textOutput("recommendation"))),
          tabsetPanel(
            tabPanel("Compare",
                     div(class = "panel",
                         plotOutput("compare_plot", height = "320px"),
                         DTOutput("compare_table"),
                         tags$hr(),
                         tags$h4("Selected-group patient-clustered GEE"),
                         div(class = "note", textOutput("regression_note")),
                         plotOutput("regression_plot", height = "300px"),
                         div(class = "note", textOutput("regression_interpretation")),
                         DTOutput("regression_table")
                     )),
            tabPanel("All Groups",
                     div(class = "panel",
                         imageOutput("all_groups_heatmap", width = "100%"),
                         DTOutput("all_table")
                     )),
            tabPanel("Evidence Notes",
                     div(class = "panel",
                         tags$p("High: n >= 100; Moderate: n >= 50; Low: n >= 30; Very low: n < 30."),
                         tags$p("Adjusted group model is used for sufficiently represented day-grade groups. Component model is used when the group-specific adjusted estimate is unavailable."),
                         tags$p("The tool supports embryo selection discussions in non-PGT single blastocyst frozen-thawed embryo transfer cycles and should be interpreted with patient-specific clinical context.")
                     ))
          ),
          div(class = "metric-grid", style = "margin-top:12px;",
              div(class = "metric", div(class = "label", "Analysis cohort"), div(class = "value", format(sum(app_group_summary$n), big.mark = ","))),
              div(class = "metric", div(class = "label", "Available groups"), div(class = "value", length(unique(app_group_summary$group)))),
              div(class = "metric", div(class = "label", "Groups n >= 50"), div(class = "value", sum(app_group_summary$n >= 50))),
              div(class = "metric", div(class = "label", "Top n>=50 LBR"), div(class = "value", paste0(max(app_group_summary$live_birth_rate[app_group_summary$n >= 50], na.rm = TRUE), "%")))
          )
        )
      )
  )
)

server <- function(input, output, session) {
  input_or_default <- function(id, default) {
    value <- input[[id]]
    if (is.null(value) || length(value) == 0 || is.na(value)) default else value
  }

  output$candidate_controls <- renderUI({
    n <- as.integer(input$candidate_n)
    lapply(seq_len(n), function(i) {
      div(class = "embryo-row",
          strong(paste("Embryo", i)),
          fluidRow(
            column(6, selectInput(paste0("day_", i), "Developmental day", choices = day_choices, selected = group_part(i, "day"))),
            column(6, selectInput(paste0("stage_", i), "Expansion stage", choices = stage_choices, selected = group_part(i, 1)))
          ),
          fluidRow(
            column(6, selectInput(paste0("icm_", i), "Inner cell mass (ICM)", choices = cell_choices, selected = group_part(i, 2))),
            column(6, selectInput(paste0("te_", i), "Trophectoderm (TE)", choices = cell_choices, selected = group_part(i, 3)))
          )
      )
    })
  })

  selected_data <- reactive({
    n <- as.integer(input$candidate_n)
    groups <- lapply(seq_len(n), function(i) {
      data.frame(
        candidate = paste("Embryo", i),
        group = find_group(
          input_or_default(paste0("day_", i), ifelse(i == 1, "D5", "D6")),
          input_or_default(paste0("stage_", i), "4"),
          input_or_default(paste0("icm_", i), "B"),
          input_or_default(paste0("te_", i), "B")
        )
      )
    }) |> bind_rows()

    groups |>
      left_join(app_group_summary, by = "group") |>
      mutate(
        rank = rank(-recommended_pred_live_birth, ties.method = "first"),
        rank_label = paste0("#", rank),
        observed_live_birth = ifelse(is.na(live_birth_rate), NA_character_, paste0(live_birth_rate, "% (", live_birth_ci_low, "-", live_birth_ci_high, ")")),
        adjusted_live_birth = ifelse(is.na(recommended_pred_live_birth), NA_character_, paste0(recommended_pred_live_birth, "%"))
      ) |>
      arrange(rank)
  })

  output$recommendation <- renderText({
    dat <- selected_data()
    validate(need(n_distinct(dat$group) == nrow(dat), "Candidate embryos must have distinct day-grade groups."))
    best <- dat |> slice_min(rank, n = 1)
    paste0(best$candidate, " (", best$group, "), adjusted estimated live-birth probability ",
           best$recommended_pred_live_birth, "%; evidence level: ", best$evidence_level, ".")
  })

  output$compare_plot <- renderPlot({
    dat <- selected_data()
    validate(need(nrow(dat) > 0, "Select embryos to compare."))
    validate(need(n_distinct(dat$group) == nrow(dat), "Candidate embryos must have distinct day-grade groups."))
    dat |>
      dplyr::select(candidate, group, recommended_pred_live_birth, clinical_rate, miscarriage_rate) |>
      tidyr::pivot_longer(c(recommended_pred_live_birth, clinical_rate, miscarriage_rate), names_to = "outcome", values_to = "rate") |>
      dplyr::mutate(
        outcome = dplyr::recode(outcome, recommended_pred_live_birth = "Live birth", clinical_rate = "Clinical pregnancy", miscarriage_rate = "Miscarriage"),
        outcome = factor(outcome, levels = c("Live birth", "Clinical pregnancy", "Miscarriage")),
        candidate_group = paste(candidate, group, sep = ": ")
      ) |>
      ggplot(aes(outcome, rate, fill = candidate_group)) +
      geom_col(position = "dodge", width = 0.72, color = "#3f4544") +
      geom_text(aes(label = paste0(rate, "%")), position = position_dodge(width = 0.72), vjust = -0.35, size = 3) +
      scale_fill_manual(values = c("#8fa39a", "#7d91a3", "#b58c87", "#9a8f9e", "#c8c1b8")) +
      scale_y_continuous(limits = c(0, 100)) +
      labs(x = NULL, y = "Rate (%)", fill = "Candidate embryo") +
      theme_classic(base_size = 12) +
      theme(panel.border = element_rect(fill = NA, color = "#3f4544"))
  })

  output$compare_table <- renderDT({
    dat <- selected_data() |>
      transmute(
        Rank = rank_label,
        Candidate = candidate,
        Group = group,
        N = n,
        `Observed live birth, % (95% CI)` = observed_live_birth,
        `Adjusted live birth` = adjusted_live_birth,
        `Clinical pregnancy, %` = clinical_rate,
        `Miscarriage, %` = miscarriage_rate,
        Evidence = evidence_level,
        Source = model_source,
        Note = evidence_note
      )
    datatable(dat, rownames = FALSE, options = list(pageLength = 5, dom = "tip", scrollX = TRUE))
  })

  selected_regression <- reactive({
    selected_groups <- unique(selected_data()$group)
    validate(need(length(selected_groups) >= 2, "Select at least two distinct embryo groups."))

    dat <- app_regression_data |>
      dplyr::filter(group %in% selected_groups) |>
      stats::na.omit()
    validate(need(nrow(dat) >= 20, "Insufficient complete cases for adjusted regression."))

    complete_counts <- dat |>
      dplyr::count(group, name = "complete_case_n") |>
      dplyr::arrange(dplyr::desc(complete_case_n), group)
    validate(need(nrow(complete_counts) >= 2, "At least two selected groups require complete-case data."))

    total_counts <- selected_data() |>
      dplyr::distinct(group, n) |>
      dplyr::filter(group %in% complete_counts$group) |>
      dplyr::arrange(dplyr::desc(n), group)
    reference <- total_counts$group[1]
    dat <- dat |>
      dplyr::mutate(
        group_model = stats::relevel(factor(group), ref = reference),
        Treatment_method_model = droplevels(Treatment_method_model),
        Embryo_derivation_model = droplevels(Embryo_derivation_model)
      ) |>
      dplyr::arrange(patient_id)

    candidate_covariates <- c(
      "Female_age", "Female_BMI", "Treatment_method_model",
      "Embryo_derivation_model", "preparation_endometrial_thickness"
    )
    usable_covariates <- candidate_covariates[vapply(
      dat[candidate_covariates],
      function(x) length(unique(x)) > 1,
      logical(1)
    )]
    model_formula <- stats::as.formula(paste(
      "live_birth ~",
      paste(c("group_model", usable_covariates), collapse = " + ")
    ))
    fit <- geepack::geeglm(
      model_formula,
      data = dat,
      id = patient_id,
      family = stats::binomial(),
      corstr = "exchangeable"
    )
    coef_table <- summary(fit)$coefficients
    group_terms <- grep("^group_model", rownames(coef_table), value = TRUE)
    result <- lapply(group_terms, function(term) {
      beta <- coef_table[term, "Estimate"]
      se <- coef_table[term, "Std.err"]
      data.frame(
        group = sub("^group_model", "", term),
        adjusted_or = exp(beta),
        conf_low = exp(beta - 1.96 * se),
        conf_high = exp(beta + 1.96 * se),
        p_value = coef_table[term, "Pr(>|W|)"]
      )
    }) |>
      dplyr::bind_rows() |>
      dplyr::filter(
        is.finite(adjusted_or), is.finite(conf_low), is.finite(conf_high), is.finite(p_value),
        adjusted_or > 0, conf_low > 0, conf_high > 0
      ) |>
      dplyr::left_join(complete_counts, by = "group") |>
      dplyr::bind_rows(data.frame(
        group = reference, adjusted_or = 1, conf_low = 1, conf_high = 1,
        p_value = NA_real_,
        complete_case_n = complete_counts$complete_case_n[complete_counts$group == reference]
      )) |>
      dplyr::mutate(reference = group == reference) |>
      dplyr::arrange(dplyr::desc(reference), dplyr::desc(complete_case_n))

    list(
      result = result,
      reference = reference,
      model_n = nrow(dat),
      model_patient_n = dplyr::n_distinct(dat$patient_id),
      events = sum(dat$live_birth),
      adjusted_terms = usable_covariates
    )
  })

  output$regression_note <- renderText({
    reg <- selected_regression()
    paste0(
      "The regression dataset contains only cycles from the currently selected embryo groups. ",
      reg$reference, " is the reference because it has the largest total cycle count among the selected groups. ",
      adjustment_label, " Covariates without variation in the selected dataset are omitted. Model cycles n=", format(reg$model_n, big.mark = ","),
      "; patients n=", format(reg$model_patient_n, big.mark = ","),
      "; live births n=", format(reg$events, big.mark = ","), "."
    )
  })

  output$regression_plot <- renderPlot({
    dat <- selected_regression()$result |>
      dplyr::filter(!reference) |>
      dplyr::mutate(group = stats::reorder(group, adjusted_or))
    validate(need(nrow(dat) > 0, "No comparison coefficient is available."))
    ggplot(dat, aes(adjusted_or, group)) +
      geom_vline(xintercept = 1, linetype = 2, color = "#3f4544") +
      geom_errorbarh(aes(xmin = conf_low, xmax = conf_high), height = 0.16, color = "#3f4544") +
      geom_point(size = 3.2, shape = 21, fill = "#9a8f9e", color = "#3f4544") +
      scale_x_log10() +
      labs(
        x = "Population-averaged adjusted odds ratio for live birth (log scale)",
        y = NULL,
        title = paste0("Reference: ", selected_regression()$reference)
      ) +
      theme_classic(base_size = 12) +
      theme(panel.border = element_rect(fill = NA, color = "#3f4544"))
  })

  output$regression_interpretation <- renderText({
    reg <- selected_regression()
    x <- reg$result |> filter(!reference)
    paste(vapply(seq_len(nrow(x)), function(i) {
      direction <- if (x$adjusted_or[i] > 1) "higher" else "lower"
      significant <- x$conf_low[i] > 1 || x$conf_high[i] < 1
      paste0("Compared with ", reg$reference, ", ", x$group[i], " had ", direction,
             " population-averaged adjusted odds of live birth (OR ", sprintf("%.2f", x$adjusted_or[i]), ", 95% CI ",
             sprintf("%.2f-%.2f", x$conf_low[i], x$conf_high[i]), "); the difference was ",
             ifelse(significant, "statistically significant.", "not statistically significant."))
    }, character(1)), collapse = " ")
  })

  output$regression_table <- renderDT({
    dat <- selected_regression()$result |>
      dplyr::transmute(
        Group = group,
        `Included cycles, n` = complete_case_n,
        `Adjusted OR` = ifelse(reference, "Reference", sprintf("%.2f", adjusted_or)),
        `95% CI` = ifelse(reference, "-", sprintf("%.2f-%.2f", conf_low, conf_high)),
        P = ifelse(reference, "-", ifelse(p_value < 0.001, "<0.001", sprintf("%.3f", p_value)))
      )
    datatable(dat, rownames = FALSE, options = list(dom = "t", scrollX = TRUE))
  })

  output$all_table <- renderDT({
    dat <- app_group_summary |>
      transmute(
        Group = group,
        Day = Embryo_day,
        Grade = Embryo_morphology_clean,
        N = n,
        `Live births` = live_birth_n,
        `Observed LBR, %` = live_birth_rate,
        `95% CI low` = live_birth_ci_low,
        `95% CI high` = live_birth_ci_high,
        `Adjusted LBR, %` = recommended_pred_live_birth,
        `Clinical pregnancy, %` = clinical_rate,
        `Miscarriage, %` = miscarriage_rate,
        Evidence = evidence_level,
        Source = model_source
      )
    datatable(dat, rownames = FALSE, filter = "top", options = list(pageLength = 15, scrollX = TRUE))
  })

  output$all_groups_heatmap <- renderImage({
    list(src = normalizePath(file.path("..", "figures", "03_live_birth_heatmap.png")),
         contentType = "image/png", alt = "Live birth heatmap for all observed blastocyst groups")
  }, deleteFile = FALSE)
}

shinyApp(ui, server)
