#!groovy​

// The 'xchem' Fragalysis Jenkinsfile.

pipeline {

  agent { label 'buildah-slave' }

  environment {
    // Slack channel for all notifications
    SLACK_BUILD_CHANNEL = 'dls-builds'
    // Slack channel to be used for errors/failures
    SLACK_ALERT_CHANNEL = 'dls-alerts'
  }

  stages {

    stage('Inspect') {
      steps {
        slackSend channel: "#${SLACK_BUILD_CHANNEL}",
                  message: "${JOB_NAME} build ${BUILD_NUMBER} - starting..."
        echo "Inspecting..."
        sh "false"
      }
    }

  }

  // Post-job actions.
  // See https://jenkins.io/doc/book/pipeline/syntax/#post
  post {

    success {
      slackSend channel: "#${SLACK_BUILD_CHANNEL}",
                color: 'good',
                message: "${JOB_NAME} build ${BUILD_NUMBER} - complete"
    }

    failure {
      slackSend channel: "#${SLACK_BUILD_CHANNEL}",
                color: 'danger',
                message: "${JOB_NAME} build ${env.BUILD_NUMBER} - failed (${BUILD_URL})"
      slackSend channel: "#${SLACK_ALERT_CHANNEL}",
                color: 'danger',
                message: "${JOB_NAME} build ${BUILD_NUMBER} - failed (${BUILD_URL})"
    }

    fixed {
      slackSend channel: "#${SLACK_ALERT_CHANNEL}",
                color: 'good',
                message: "${JOB_NAME} build - fixed"
    }

  }

}
