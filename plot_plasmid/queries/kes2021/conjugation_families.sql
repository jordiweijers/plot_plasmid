-- MCL families whose most common protein symbol is tra followed by one capital letter (traA, traB, ...).
-- Only proteins with a symbol are counted, not pseudogenes; a tie goes to the alphabetically first symbol.
-- REGEXP BINARY makes the match case-sensitive.
-- Queries the whole database, so the category saves its IDs at its path.
WITH symbol_counts AS (
    SELECT m.clust, p.symbol, COUNT(*) AS symbol_n
    FROM proteins p
    JOIN mcl30 m ON m.id = p.id
    WHERE p.symbol IS NOT NULL AND p.symbol != ''
    GROUP BY m.clust, p.symbol
),
ranked_symbols AS (
    SELECT clust, symbol,
           ROW_NUMBER() OVER (PARTITION BY clust ORDER BY symbol_n DESC, symbol ASC) AS rn
    FROM symbol_counts
)
SELECT clust
FROM ranked_symbols
WHERE rn = 1 AND symbol REGEXP BINARY '^tra[A-Z]$'
