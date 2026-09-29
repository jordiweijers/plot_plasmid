-- MCL families whose most common protein name contains "replication initiat" (REGEXP ignores case).
-- Only proteins with a name are counted: CDSs with 'pseudo' in their attributes are skipped.
-- A tie goes to the alphabetically first name.
-- Queries the whole database, so the category saves its IDs at its path.
WITH name_counts AS (
    SELECT m.mcl_id, f.name, COUNT(*) AS name_n
    FROM features_cds f
    JOIN mcls m ON m.cds_id = f.cds_id
    WHERE f.name IS NOT NULL AND f.name != ''
      AND COALESCE(FIND_IN_SET('pseudo', REPLACE(f.attributes, ';', ',')), 0) = 0
    GROUP BY m.mcl_id, f.name
),
ranked_names AS (
    SELECT mcl_id, name,
           ROW_NUMBER() OVER (PARTITION BY mcl_id ORDER BY name_n DESC, name ASC) AS rn
    FROM name_counts
)
SELECT mcl_id
FROM ranked_names
WHERE rn = 1 AND name REGEXP 'replication initiat'
