-- MCL families whose most common protein name contains "replication initiat" (REGEXP ignores case).
-- Only proteins with a name are counted, not pseudogenes; a tie goes to the alphabetically first name.
-- Queries the whole database, so the category saves its IDs at its path.
WITH name_counts AS (
    SELECT m.clust, p.name, COUNT(*) AS name_n
    FROM proteins p
    JOIN mcl30 m ON m.id = p.id
    WHERE p.name IS NOT NULL AND p.name != ''
    GROUP BY m.clust, p.name
),
ranked_names AS (
    SELECT clust, name,
           ROW_NUMBER() OVER (PARTITION BY clust ORDER BY name_n DESC, name ASC) AS rn
    FROM name_counts
)
SELECT clust
FROM ranked_names
WHERE rn = 1 AND name REGEXP 'replication initiat'
