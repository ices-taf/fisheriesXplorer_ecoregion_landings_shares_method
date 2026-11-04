library(icesTAF)
taf.skeleton()


library(dplyr)

# 1. Basic size check
dim(cl_old)
dim(cl_new)

# 2. Check total scientific weight is unchanged
sum(cl_old$CLscientificWeight, na.rm = TRUE)
sum(cl_new$CLscientificWeight, na.rm = TRUE)

# 3. Compare how many FAO codes are missing
sum(is.na(cl_old$CLspeciesFaoCode) | cl_old$CLspeciesFaoCode == "")
sum(is.na(cl_new$CLspeciesFaoCode) | cl_new$CLspeciesFaoCode == "")

# 4. Compare the actual FAO-code assignments
comparison <- cl_old %>%
  transmute(
    row_id = row_number(),
    old_fao = CLspeciesFaoCode
  ) %>%
  left_join(
    cl_new %>%
      transmute(
        row_id = row_number(),
        new_fao = CLspeciesFaoCode
      ),
    by = "row_id"
  ) %>%
  mutate(
    old_fao = na_if(old_fao, ""),
    new_fao = na_if(new_fao, ""),
    same = old_fao == new_fao |
      (is.na(old_fao) & is.na(new_fao))
  )

# 5. How many rows changed?
comparison %>%
  count(same)

# 6. Show only differences
comparison %>%
  filter(!same | is.na(same))


comparison <- comparison %>%
  mutate(
    status = case_when(
      is.na(old_fao) & is.na(new_fao) ~ "both missing",
      is.na(old_fao) & !is.na(new_fao) ~ "filled by new method",
      !is.na(old_fao) & is.na(new_fao) ~ "lost in new method",
      old_fao == new_fao ~ "same",
      TRUE ~ "different code"
    )
  )

comparison %>%
  count(status, sort = TRUE)


# Cases where the new method changed an existing non-missing code
comparison %>%
  filter(status == "different code") %>%
  count(old_fao, new_fao, sort = TRUE)

# Cases where the old method had a code but the new method lost it
comparison %>%
  filter(status == "lost in new method") %>%
  count(old_fao, new_fao, sort = TRUE)

comparison %>%
  filter(status == "different code")


comparison2 <- cl_old %>%
  transmute(
    row_id = row_number(),
    CLspeciesCode,
    old_fao = na_if(CLspeciesFaoCode, "")
  ) %>%
  left_join(
    cl_new %>%
      transmute(
        row_id = row_number(),
        new_fao = na_if(CLspeciesFaoCode, "")
      ),
    by = "row_id"
  ) %>%
  mutate(
    status = case_when(
      is.na(old_fao) & is.na(new_fao) ~ "both missing",
      is.na(old_fao) & !is.na(new_fao) ~ "filled by new method",
      !is.na(old_fao) & is.na(new_fao) ~ "lost in new method",
      old_fao == new_fao ~ "same",
      TRUE ~ "different code"
    )
  )

comparison2 %>%
  filter(status %in% c("lost in new method", "different code")) %>%
  count(CLspeciesCode, old_fao, new_fao, status, sort = TRUE)
