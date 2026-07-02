system(
  sprintf(
    'for inp in "%s"/preprocess_*.rds; do
      Rscript scripts/cluster-sobj.R --input "$inp" --elbow-n 20
    done
  ',
    CURRENT_OBJECT_DIR
  )
)
