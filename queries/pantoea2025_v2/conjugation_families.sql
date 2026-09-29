-- MCL families whose most common protein symbol is tra followed by one capital letter (traA, traB, ...).
-- Only proteins with a symbol are counted: CDSs with 'pseudo' in their attributes are skipped.
-- A tie goes to the alphabetically first symbol. REGEXP BINARY makes the match case-sensitive.
-- Queries the whole database, so the category saves its IDs at its path.
WITH symbol_counts AS (
    SELECT m.mcl_id, f.symbol, COUNT(*) AS symbol_n
    FROM features_cds f
    JOIN mcls m ON m.cds_id = f.cds_id
    WHERE f.symbol IS NOT NULL AND f.symbol != ''
      AND COALESCE(FIND_IN_SET('pseudo', REPLACE(f.attributes, ';', ',')), 0) = 0
    GROUP BY m.mcl_id, f.symbol
),
ranked_symbols AS (
    SELECT mcl_id, symbol,
           ROW_NUMBER() OVER (PARTITION BY mcl_id ORDER BY symbol_n DESC, symbol ASC) AS rn
    FROM symbol_counts
)
SELECT mcl_id
FROM ranked_symbols
WHERE rn = 1 AND symbol REGEXP BINARY '^tra[A-Z]$'
