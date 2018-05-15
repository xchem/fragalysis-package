docker run --publish=7474:7474 --publish=7687:7687  --volume=$PWD:/data --env=NEO4J_AUTH=none   neo4j
