Shiny.addCustomMessageHandler("testmessageHandlerForJS",
  function(message) {
    alert(JSON.stringify(message));
  }
);
